package heim

import heim.tasks.DataFetchTaskFactory
import heim.tasks.TaskState
import org.transmartproject.db.dataquery.highdim.mrna.MrnaTestData

import static org.hamcrest.MatcherAssert.assertThat
import static org.hamcrest.Matchers.*

/**
 * Created by piotrzakrzewski on 13/10/15.
 */
class ScriptExecutionControllerTests extends BaseAPITestCase {


    private static final String SAMPLE_DATA_TO_PASS = 'sample data to pass'
    private static final String SAMPLE_DATA_TO_WRITE = 'sample data to write'
    private static final String DATA_FETCH_TASK_TYPE = DataFetchTaskFactory.FETCH_DATA_TASK_NAME


    private Map runRScript(String name, Map arguments = [:], String workflow = '_func_test') {
        String sessionId = createSession(workflow)
        post '/ScriptExecution/run', {
            body json: [
                    sessionId: sessionId,
                    taskType:  name,
                    arguments: arguments
            ]
        }

        assertStatus 200
        assertThat JSON, hasEntry(
                equalTo('executionId'),
                isA(String))

        String taskId = JSON.executionId
        get '/ScriptExecution/status?' + buildQueryParameters(
                sessionId:         sessionId,
                executionId:       taskId,
                waitForCompletion: true,
        )

        [sessionId: sessionId, taskId: taskId]
    }

    private void fetchTestMRNA(sessionId){
        post '/ScriptExecution/run', {
            body json: [
                    sessionId: sessionId,
                    taskType: DATA_FETCH_TASK_TYPE,
                    arguments: [
                            conceptKeys: [test_label: '\\\\i2b2 main\\foo\\study1\\bar\\'],
                            dataType: 'mrna',
                            resultInstanceIds: [], //TODO: the test data has no result set for the mrna patients
                            assayConstraints: [
                                    trial_name: [name: MrnaTestData.TRIAL_NAME],
                            ],
                            dataConstraints: [
                                    genes: [ names: ['BOGUSRQCD1'] ]
                            ],
                            projection: 'zscore',
                    ]
            ]
        }
        String taskId = JSON.executionId
        get '/ScriptExecution/status?' + buildQueryParameters(
                sessionId:         sessionId,
                executionId:       taskId,
                waitForCompletion: true,
        )

        assertStatus(200)
    }

    private Map runHeatmap() {
        String sessionId = createSession('heatmap')
        fetchTestMRNA(sessionId)
        post '/ScriptExecution/run', {
            body json: [
                    sessionId: sessionId,
                    taskType: 'run',
                    arguments: []
            ]
        }

        assertStatus 200
        assertThat JSON, hasEntry(is('executionId'), isA(String))

        String taskId = JSON.executionId
        get '/ScriptExecution/status?' + buildQueryParameters(
                sessionId:         sessionId,
                executionId:       taskId,
                waitForCompletion: true,
        )
    }

    void testHeatmap() {
        return
        runHeatmap()

        assertStatus 200

        assertThat JSON, allOf(
                hasEntry(is('state'), is(TaskState.FINISHED.toString())),
                hasEntry(is('result'), allOf(
                        hasEntry(is('successful'), is(true)),
                        hasEntry(is('exception'), is(null)), // not is(nullValue()); it's a JSONObject.Null
                        hasEntry(is('artifacts'), isA(Map))
                )))

        Map artifacts = JSON.result.artifacts
        assertThat artifacts['files'], hasItem('heatmap.json')
    }

    void testTouch() {
        String sessionId = createSession('_func_test')
        post '/RSession/touch', {
            body json: [
                    sessionId: sessionId,
            ]
        }

        assertStatus 204
    }


    void testRun() {
        def args = [
                data_to_pass: SAMPLE_DATA_TO_PASS,
                data_to_write: SAMPLE_DATA_TO_WRITE,
        ]
        runRScript('sample', args)

        assertStatus 200
        assertThat JSON, allOf(
                hasEntry(is('state'), is(TaskState.FINISHED.toString())),
                hasEntry(is('result'), allOf(
                        hasEntry(is('successful'), is(true)),
                        hasEntry(is('exception'), is(null)), // not is(nullValue()); it's a JSONObject.Null
                        hasEntry(is('artifacts'), isA(Map))
                )))

        Map artifacts = JSON.result.artifacts
        assertThat artifacts['files'], contains('test_file')
        assert artifacts['a'] == (1 as Double)
        assert artifacts['b'] == 'foobar'
        assert artifacts['passed'] == SAMPLE_DATA_TO_PASS
    }

    void testSourcing() {
        runRScript('sourcing')
        assertStatus 200
        Map artifacts = JSON.result.artifacts
        assert artifacts['shouldBeTest'] == 'test'
    }

    void testCoreUtils() {
        runRScript('autosourcing')
        assertStatus 200
    }

    void testRetrieveFile() {
        def args = [
                data_to_pass: SAMPLE_DATA_TO_PASS,
                data_to_write: SAMPLE_DATA_TO_WRITE,
        ]
        def m = runRScript('sample', args)

        get '/ScriptExecution/downloadFile?' + buildQueryParameters(
                sessionId:         m.sessionId,
                executionId:       m.taskId,
                filename:          JSON.result.artifacts.files[0],
        )

        assertThat content, containsString(SAMPLE_DATA_TO_WRITE)
    }

    void testRetrieveFileBadArguments() {
        def args = [
                data_to_pass: SAMPLE_DATA_TO_PASS,
                data_to_write: SAMPLE_DATA_TO_WRITE,
        ]
        runRScript('sample',args)

        // don't give any arguments
        get '/ScriptExecution/downloadFile'

        assertStatus 400
    }

}
