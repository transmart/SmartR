package heim.session

import grails.util.Environment
import groovy.util.logging.Log4j
import heim.SmartRExecutorService
import heim.SmartRRuntimeConstants
import heim.jobs.JobInstance
import heim.tasks.JobTasksService
import heim.tasks.TaskResult
import heim.tasks.TaskState
import org.springframework.beans.factory.DisposableBean
import org.springframework.beans.factory.annotation.Autowired
import org.springframework.stereotype.Component
import org.transmartproject.core.exceptions.InvalidArgumentsException
import org.transmartproject.core.exceptions.NoSuchResourceException
import org.transmartproject.core.users.User

import javax.annotation.PostConstruct
import java.util.concurrent.Callable
import java.util.concurrent.atomic.AtomicReference
import java.util.concurrent.locks.Lock
import java.util.concurrent.locks.ReadWriteLock
import java.util.concurrent.locks.ReentrantReadWriteLock

/**
 * Public service with the controller is to communicate.
 */
@Component
@Log4j
class SessionService implements DisposableBean {

    private static final String[] HIDDEN_WORKFLOWS = []

    private static final int COLLECTION_INTERVAL = 5 * 60 * 1000 // 5 min
    private static final int SESSION_LIFESPAN = 10 * 60 * 1000 // 10 min

    private final ReadWriteLock sessionBookKeepingLock = new ReentrantReadWriteLock()
    // synchronizes the next two variables:
    private final Map<UUID, SessionContext> currentSessions = new HashMap<>()
    private final Set<UUID> sessionsShuttingDown = [] as Set

    private final AtomicReference<Timer> gcTimer = new AtomicReference()

    @Autowired
    SmartRRuntimeConstants constants

    @Autowired
    SmartRExecutorService smartRExecutorService

    @Autowired
    private JobTasksService jobTasksService // session scoped

    @Autowired
    private JobInstance jobInstance // session scoped

    @Autowired
    private SessionFiles sessionFiles

    @PostConstruct
    private startCollecting() {
        gcTimer.set(
                new Timer('SmartR-Session-GC', true).schedule(
                        this.&garbageCollection as TimerTask,
                        COLLECTION_INTERVAL, COLLECTION_INTERVAL))
    }

    private <T> T withSessionBookKeepingLock(boolean write, Closure<T> closure) {
        Lock lock
        if (write) {
            lock = sessionBookKeepingLock.writeLock()
        } else {
            lock = sessionBookKeepingLock.readLock()
        }
        lock.lock()
        try {
            closure.call()
        } finally {
            lock.unlock()
        }
    }

    List<String> availableWorkflows() {
        File dir = constants.pluginScriptDirectory
        dir.listFiles({ File f ->
            f.isDirectory() && !f.name.startsWith('_') &&
                    (Environment.current == Environment.DEVELOPMENT || !(f.name in HIDDEN_WORKFLOWS))
        } as FileFilter)*.name.sort()
    }

    UUID createSession(User user, String workflowType) {
        SessionContext newSession = new SessionContext(user, workflowType)
        log.debug("Created session with id ${newSession.sessionId} and " +
                "workflow type $workflowType")

        withSessionBookKeepingLock(true) {
            currentSessions[newSession.sessionId] = newSession
        }

        newSession.sessionId
    }

    void destroySession(UUID sessionId) {
        if (sessionId == null) {
            throw new NullPointerException("sessionId must be given")
        }

        SessionContext sessionContext
        withSessionBookKeepingLock(true) {
            sessionContext = fetchOperationalSessionContext(sessionId)
            if (!sessionContext) {
                throw new InvalidArgumentsException(
                        "No such operational session: $sessionId")
            }

            sessionsShuttingDown << sessionId
        }

        smartRExecutorService.submit({
            log.debug("Started callable for destroying session $sessionId")
            SmartRSessionSpringScope.withActiveSession(sessionContext) {
                log.debug(
                        "Running destruction callbacks for session $sessionId")
                sessionContext.destroy()
                log.debug("Finished running destruction " +
                        "callbacks for session $sessionId")

                withSessionBookKeepingLock(true) {
                    currentSessions.remove(sessionId)
                    sessionsShuttingDown.remove(sessionId)
                }
            }
        } as Callable<Void>)

        log.debug("Submmitted session $sessionId for destruction")
    }

    UUID createTask(Map<String, Object> arguments,
                    UUID sessionId,
                    String taskType) {
        doWithSession(sessionId) {
            def task = jobInstance.createTask(taskType, arguments)
            jobTasksService.submitTask(task)
            task.uuid
        }

    }

    Map<String, Object> getTaskData(UUID sessionUUID,
                                    UUID taskUUID,
                                    boolean waitForCompletion = false) {
        doWithSession(sessionUUID) {
            def state = jobTasksService.getTaskState(taskUUID)
            if (!state) {
                throw new NoSuchResourceException(
                        "No task $taskUUID for session $sessionUUID")
            }

            TaskResult result
            if (waitForCompletion) {
                result = jobTasksService.getTaskResultFuture(taskUUID).get()
                if (state == TaskState.QUEUED || state == TaskState.RUNNING) {
                    state = jobTasksService.getTaskState(taskUUID)
                }
            } else {
                // not great code, this may be out of sync with the state we got before
                result = jobTasksService.getTaskResult(taskUUID)
            }

            [
                    state : state,
                    result: result, // null if the task has not finished
            ]
        }
    }

    public <T> T doWithSession(UUID sessionUUID, Closure<T> closure) {
        SessionContext context = fetchOperationalSessionContext(sessionUUID)
        if (context == null) {
            throw new NoSuchResourceException(
                    "No such operational session: $sessionUUID")
        }
        context.updateLastModified()
        SmartRSessionSpringScope.withActiveSession(
                context, closure)
    }

    File getFile(UUID sessionUUID, UUID taskId, String filename) {
        def res = doWithSession(sessionUUID) {
            sessionFiles.get(taskId, filename)
        }
        if (res == null) {
            throw new NoSuchResourceException(
                    "No file '$filename' for sesison $sessionUUID and task $taskId")
        }
        res
    }

    void removeAllFiles(UUID sessionUUID) {
        doWithSession(sessionUUID) {
            sessionFiles.removeAll()
        }
    }

    boolean isSessionActive(UUID sessionId) {
        withSessionBookKeepingLock(false) {
            currentSessions.containsKey(sessionId) &&
                    !(sessionId in sessionsShuttingDown)
        }
    }

    void touchSession(sessionId) {
        try {
            doWithSession(sessionId) {}
        }
        catch (NoSuchResourceException e) {
            log.warn("Attempted to touch non-existent or shutting down " +
                    "session: $sessionId. This is normal if the session was " +
                    "destroyed with tasks running.")
        }
    }

    private SessionContext fetchOperationalSessionContext(UUID sessionId) {
        withSessionBookKeepingLock(false) {
            SessionContext sessionContext = fetchSessionContext(sessionId)
            if (sessionContext == null) {
                return null
            }
            if (sessionId in sessionsShuttingDown) {
                log.warn("Session is shutting down: $sessionId")
                return null
            }
            sessionContext
        }
    }

    // no locking, but only called from fetchOperationalSessionContext
    private SessionContext fetchSessionContext(UUID sessionId) {
        SessionContext sessionContext = currentSessions[sessionId]
        if (!sessionContext) {
            log.warn("No such session: $sessionId")
            return null
        }
        sessionContext
    }


    void garbageCollection() {
        log.debug('Started session garbage collecting')

        // because we're holding the lock, no actual destruction will start
        // until this method returns (this function only schedules destruction
        // in other threads)
        withSessionBookKeepingLock(true) {
            currentSessions.each {
                sessionId, sessionContext ->
                    def lastActivity = sessionContext.lastActive
                    if (isStale(sessionId, lastActivity)) {
                        if (!sessionsShuttingDown.contains(sessionId)) {
                            destroySession(sessionId)
                            log.info("Terminated session: ${sessionId} " +
                                    "due to inactivity.")
                        } else {
                            log.info("Session $sessionId is stale, but " +
                                    "it's already shutting down")
                        }
                    }
            }
        }
        log.debug('Finished session garbage collecting')
    }

    private boolean isStale(UUID sessionId, Date lastTouched) {
        SessionContext context = fetchOperationalSessionContext(sessionId)
        if (context == null) {
            log.debug("Session $sessionId shut down between call to " +
                    "getCurrentSessions() and fetchOperationalSessionContext()?")
            return false
        }

        SmartRSessionSpringScope.withActiveSession(context) {
            Date now = new Date()

            def delta = now.time - lastTouched.time
            boolean deltaPassed = delta > SESSION_LIFESPAN
            boolean hasActiveTasks = jobTasksService.hasActiveTasks()

            boolean res = deltaPassed && !hasActiveTasks
            if (res && log.debugEnabled) {
                log.debug("Session $sessionId deemed stale " +
                        "($delta > $SESSION_LIFESPAN and no active tasks)")
            } else if (!res && log.debugEnabled) {
                log.debug("Session $sessionId deemed active " +
                        "($delta <= $SESSION_LIFESPAN) or (active tasks: $hasActiveTasks)")
            }

            res
        }
    }

    @Override
    void destroy() throws Exception {
        gcTimer.get().cancel()
    }
}
