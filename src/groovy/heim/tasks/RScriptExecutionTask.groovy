package heim.tasks

import com.google.common.collect.ImmutableList
import com.google.common.collect.ImmutableMap
import groovy.transform.TypeChecked
import groovy.util.logging.Log4j
import heim.SmartRRuntimeConstants
import heim.rserve.RFunctionArg
import heim.rserve.RScriptOutputManager
import heim.rserve.RServeSession
import heim.rserve.RUtil
import heim.session.SessionFiles
import org.rosuda.REngine.REXP
import org.rosuda.REngine.Rserve.RConnection

import static heim.rserve.RUtil.runRCommand

@Log4j
@TypeChecked
// TODO: convert to prototype bean
class RScriptExecutionTask extends AbstractTask {

    UUID sessionId
    RServeSession rServeSession
    SessionFiles sessionFiles
    File fileToLoad // absolute
    String function = 'main'
    List<RFunctionArg> arguments
    SmartRRuntimeConstants constants

    @Override
    TaskResult call() throws Exception {
        rServeSession.doWithRConnection { RConnection conn ->
            def outputManager = new RScriptOutputManager(conn, sessionId, uuid, constants)

            REXP res = callR(conn)
            List<File> fileNames = outputManager.downloadFiles()

            convertSuccessfulResult(res, fileNames)
        }
    }

    private void injectScriptDir(RConnection conn) {
        def remoteScriptDir = fileToLoad.getParentFile().getParentFile()  // We need the HeimScripts root
        String path = RUtil.escapeRStringContent(remoteScriptDir.absolutePath)
        runRCommand(conn, "remoteScriptDir <- \"$path\"")
    }

    private void sourceCoreUtils(RConnection conn) {
        def remoteScriptDirRoot = fileToLoad.getParentFile().getParentFile()
        File remoteCoreUtilsDir = new File(remoteScriptDirRoot, '_core')
        File coreUtilsIndex = new File(remoteCoreUtilsDir, 'index.R')
        runRCommand(conn, "source('" +
                "${RUtil.escapeRStringContent(coreUtilsIndex.absolutePath)}')")
    }

    private REXP callR(RConnection conn) {
        runRCommand(conn, 'library(jsonlite)')
        injectScriptDir(conn)
        sourceCoreUtils(conn)
        runRCommand(conn, "source('" +
                "${RUtil.escapeRStringContent(fileToLoad.toString())}')")
        def namedArguments = arguments.collect { RFunctionArg arg ->
            "${arg.name}=${arg.asRExpression()}"
        }
        runRCommand(conn, "main(${namedArguments.join(', ')})")
    }

    private TaskResult convertSuccessfulResult(REXP result,
                                               List<File> files) {
        def builder = ImmutableMap.builder()
        builder.putAll(massage(result))
        builder.put('files', ImmutableList.copyOf(files*.name))

        files.each { File f ->
            sessionFiles.add(uuid, f.name, f)
        }
        new TaskResult(
                successful: true,
                artifacts: builder.build(),
        )
    }

    private ImmutableMap<String, Object> massage(REXP result) {
        def asNative = result.asNativeJavaObject()
        if (!(asNative instanceof Map)) {
            asNative = [(asNative): 'value']
        }

        /* list(a=1, b='a') results in a map where the keys (!) are arrays
         * with one element with value 1.0 and 'a' and the values are the names
         * of the list. Invert this and unwrap the values
         */
        def builder = ImmutableMap.builder()
        asNative.each { k, v ->
            if (k == []) {
                return
            }
            assert k.getClass().isArray() && ((Object[])k).length == 1
            assert v instanceof String
            builder.put(v,  ((Object[])k)[0])
        }

        builder.build()
    }

    @Override
    void close() throws Exception {
        // nothing to do
    }
}
