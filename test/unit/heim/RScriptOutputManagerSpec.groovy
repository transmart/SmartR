package heim

import grails.test.mixin.TestMixin
import grails.test.mixin.support.GrailsUnitTestMixin
import heim.rserve.RScriptOutputManager
import org.rosuda.REngine.Rserve.RConnection
import spock.lang.Specification

/**
 * See the API for {@link grails.test.mixin.support.GrailsUnitTestMixin} for usage instructions
 */
@TestMixin(GrailsUnitTestMixin)
class RScriptOutputManagerSpec extends Specification {

    def setup() {

    }

    def cleanup() {
        new File(SmartRRuntimeConstants.instance.baseDir, "bogusfolder/init/plot.png").delete() // Access to a property of SmartRRuntimeConstants and "instance" mechanism inspired by http://stackoverflow.com/a/13846141/535203
    }

    void "test getScriptOutput() returning single plot"() {
        given:
        def conn = new RConnection()
        conn.eval('png("plot.png")')
        conn.eval('plot(c(0,1),c(1,1))')
        conn.eval('dev.off()')
        def sessionId = UUID.randomUUID()
        def taskId = UUID.randomUUID()
        def mngr = new RScriptOutputManager(conn, sessionId, taskId, SmartRRuntimeConstants.instance)

        when: ""
            def result = mngr.downloadFiles()
        then: "Result is a list containing plot.png and only that one file."
            result
            result.size() == 1
            result[0] instanceof File
            result[0].getAbsolutePath() == SmartRRuntimeConstants.instance.baseDir.toString()+"/$sessionId/$taskId/plot.png"
    }
}
