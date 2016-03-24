package smartR.plugin

import grails.converters.JSON
import heim.session.SessionService
import org.codehaus.groovy.grails.web.json.JSONArray
import org.codehaus.groovy.grails.web.json.JSONObject

class SmartRController {

    SessionService sessionService

    static layout = 'smartR'

    def index() {
        [ scriptList: sessionService.availableWorkflows()]
    }

    /**
    *   Called to get the path to smartR.js such that the plugin can be loaded in the datasetExplorer
    */
    def loadScripts = {

        // list of required javascript files
        def scripts = [servletContext.contextPath + pluginContextPath + '/js/smartR/smartR.js']

        // list of required css files
        def styles = []

        JSONObject result = new JSONObject()
        JSONArray rows = new JSONArray()

        // for all js files
        for (file in scripts) {
            def m = [:]
            m["path"] = file.toString()
            m["type"] = "script"
            rows.put(m);
        }

        // for all css files
        for (file in styles) {
            def n = [:]
            n["path"] = file.toString()
            n["type"] = "css"
            rows.put(n);
        }

        result.put("success", true)
        result.put("totalCount", scripts.size())
        result.put("files", rows)

        render result as JSON;
    }

    /**
     * Get smart-r plugin context path
     */
    def smartRContextPath = {
        render servletContext.contextPath + pluginContextPath as String;
    }

    def createPDMapLayout = {
        
    }
}
