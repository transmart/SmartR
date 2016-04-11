package smartR.plugin

import grails.converters.JSON
import heim.session.SessionService
import org.codehaus.groovy.grails.web.json.JSONArray
import org.codehaus.groovy.grails.web.json.JSONObject
import groovyx.net.http.HTTPBuilder
import static groovyx.net.http.Method.POST
import static groovyx.net.http.ContentType.*

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

    /**
     * Send gene expression to MINERVA server to generate a custom layout
     *
     * @param params map which should contain
     *     `identifier`:       id of custom MINERVA layout
     *     `login`:            login name
     *     `password`:         password
     *     `model`:            the name of the MINERVA model, eg. 'pdmap_dec15'
     *     `expression_value`: tab-separated values (as sting) containing `gene name` and `expression value`
     * @return json document containing the key `success`, and if its value is true the keys `url` and `rawData` (html response from MINERVA)
     */
    def createMinervaLayout = {
        def url = grailsApplication.config.minervaConnector.url ?: 'http://localhost:8080'
        def path = grailsApplication.config.minervaConnector.createCustomLayoutPath ?: '/minerva/galaxy.xhtml'
        def http = new HTTPBuilder(url)
        JSONObject result = new JSONObject()

        http.request( POST ) {
            uri.path = path
            contentType = TEXT
            send URLENC, params

            response.success = { resp, data ->
                log.debug "POST success (" + url + path + "), response status: ${resp.statusLine} " + resp.statusLine.statusCode
                result.put("success", true)
                result.put("url", url)
                result.put("rawData", data.getText())
            }

            response.failure = { resp ->
                log.warn "POST failed (" + url + path + "), response status: ${resp.statusLine} " + resp.statusLine.statusCode
                result.put("success", false)
                result.put("status", resp.status)
            }
        }

        render result as JSON;
    }
}
