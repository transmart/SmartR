// Copyright (C) 2016 ITTM S.A., Nils Christian <nils.christian@ittm-solutions.com>

package smartR.plugin

import grails.converters.JSON
import com.google.gson.Gson
import grails.util.Holders
import com.ittm_solutions.ipacore.IpaApiCached
import com.ittm_solutions.ipacore.IpaApiException
import com.ittm_solutions.ipacore.IpaApiAuthenticationFailedException

class IpaConnectorController {
    static scope = "singleton"

    def projectName = 'transmartipaapi'

    def userToIpaApiCached = [:]

    ////////////////////////////////////////////////////////////////
    // private

    /**
     * Get IpaApiCached instance.
     **/
    private def synchronized getIpaApiCached(username, password) {
        // FIXME this should create one instance per
        // (username,password,ipaApiUrl,projectName) tuple, not per username
        // (the password should be hashed)

        // NOTE this would be an alternative if one tranSMART instance
        // uses a shared IPA login (check whether your IPA license
        // allows this)
        // username = grailsApplication.config.ipaConnector.username
        // password = grailsApplication.config.ipaConnector.password

        // in case we are being handed `null` instead of string (user
        // didn't fill in fields) just use empty strings
        if (!username) {
            username = ''
        }
        if (!password) {
            password = ''
        }
        if (! userToIpaApiCached.containsKey(username)) {
            log.debug "Creating new IpaApiCached for " + username
            userToIpaApiCached[username] = new IpaApiCached(username,password,projectName)
        }
        return userToIpaApiCached[username]
    }

    private def synchronized removeIpaApiCached(username) {
        if (userToIpaApiCached.containsKey(username)) {
            userToIpaApiCached.remove(username)
        }
    }

    ////////////////////////////////////////////////////////////////

    /**
    *  Renders the default view.
    */
    def index = {
      [model: ['nada']]
    }

    /**
     * Uploads a dataset and initiates an analysis on the server.
     * <p>
     * See IpaApiCached.dataAnalysis for details.
     */
    def dataAnalysis() {
        def analysisParams = request.JSON
        def myipaapicached = getIpaApiCached(analysisParams.username, analysisParams.password)
        try {
            def retval = [
                startedNewCalculation: myipaapicached.dataAnalysis(analysisParams.analysisName,
                                                                   analysisParams.geneIdType,analysisParams.geneId,
                                                                   analysisParams.expValueType,analysisParams.expValue,
                                                                   analysisParams.expValueType2,analysisParams.expValue2,
                                                                   analysisParams.applicationName),
            ]
            render retval as JSON
        } catch (IpaApiAuthenticationFailedException e) {
            removeIpaApiCached(analysisParams.username)
            render(status: 401, text: 'authentication failed')
        } catch (IpaApiException|Exception e) {
            removeIpaApiCached(analysisParams.username)
            log.error("IpaApiCached exception:",e)
            render(status: 500, text: e.getMessage())
        }
    }

    /**
     * Returns the id of the analysis with the given name if it exists.
     * <p>
     * See IpaApiCached.analysisId for details.
     */
    def analysisId() {
        try {
            def myipaapicached = getIpaApiCached(request.JSON.username, request.JSON.password)
            def retval = [
                analysisId: myipaapicached.analysisId(request.JSON.analysisName),
            ]
            render retval as JSON
        } catch (IpaApiAuthenticationFailedException e) {
            removeIpaApiCached(analysisParams.username)
            render(status: 401, text: 'authentication failed')
        } catch (IpaApiException|Exception e) {
            removeIpaApiCached(analysisParams.username)
            log.error("IpaApiCached exception:",e)
            render(status: 500, text: e.getMessage())
        }
    }

    /**
     * Exports and returns the given analysis.
     * <p>
     * See IpaApiCached.exportIngenuityResults for details.
     */
    def exportIngenuityResults() {
        try {
            def myipaapicached = getIpaApiCached(request.JSON.username, request.JSON.password)
            def ipaResults = myipaapicached.exportAnalysis(request.JSON.analysisId)
            Gson gson = new Gson();
            def json = gson.toJson(ipaResults);
            render contentType: 'application/json', text: json
        } catch (IpaApiAuthenticationFailedException e) {
            removeIpaApiCached(analysisParams.username)
            render(status: 401, text: 'authentication failed')
        } catch (IpaApiException|Exception e) {
            removeIpaApiCached(analysisParams.username)
            log.error("IpaApiCached exception:",e)
            render(status: 500, text: e.getMessage())
        }
    }
}
