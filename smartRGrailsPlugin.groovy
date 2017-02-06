import grails.util.Environment
import heim.rserve.RScriptsSynchronizer
import org.codehaus.groovy.grails.plugins.GrailsPluginUtils
import org.springframework.stereotype.Component
import heim.SmartRRuntimeConstants

class smartRGrailsPlugin {

    public static final String DEFAULT_REMOTE_RSCRIPTS_DIRECTORY = '/tmp/smart_r_scripts'
    public static final String TRANSMART_EXTENSIONS_REGISTRY_BEAN_NAME = 'transmartExtensionsRegistry'

    // the plugin version
    def version = "1.0-STABLE-SNAPSHOT"
    // the version or versions of Grails the plugin is designed for
    def grailsVersion = "2.3 > *"
    // resources that are excluded from plugin packaging
    def pluginExcludes = [
        "grails-app/views/error.gsp"
    ]

    // TODO Fill in these fields
    def title = "SmartR Plugin" // Headline display name of the plugin
    def author = ""
    def authorEmail = ""
    def description =
            '''
            Brief summary/description of the plugin.
            '''

    // URL to the plugin's documentation
    def documentation = ""

    // Extra (optional) plugin metadata

    // License: one of 'APACHE', 'GPL2', 'GPL3'
//    def license = "APACHE"

    // Details of company behind the plugin (if there is one)
//    def organization = [ name: "My Company", url: "http://www.my-company.com/" ]


    // Any additional developers beyond the author specified above.
//    def developers = [ [ name: "Joe Bloggs", email: "joe@bloggs.net" ]]

    // Location of the plugin's issue tracker.
//    def issueManagement = [ system: "JIRA", url: "http://jira.grails.org/browse/GPMYPLUGIN" ]

    // Online location of the plugin's browseable source code.
//    def scm = [ url: "http://svn.codehaus.org/grails-plugins/" ]

    def doWithWebDescriptor = { xml ->
        // TODO Implement additions to web.xml (optional), this event occurs before
    }

    def doWithSpring = {
        xmlns context:"http://www.springframework.org/schema/context"

        context.'component-scan'('base-package': 'heim') {
            context.'include-filter'(
                    type:       'annotation',
                    expression: Component.canonicalName)
        }
    }

    def doWithApplicationContext = { ctx ->
        def config = application.config
        SmartRRuntimeConstants constants = ctx.getBean(SmartRRuntimeConstants)

        File smartRDir = GrailsPluginUtils.getPluginDirForName('smart-r')?.file
        if (!smartRDir) {
            String pluginPath = ctx.pluginManager.allPlugins.find {
                it.name == 'smartR'
            }.pluginPath

            smartRDir = ctx.getResource(pluginPath).file
        } else {
            smartRDir = new File(smartRDir, 'web-app')
        }
        if (!smartRDir) {
            throw new RuntimeException('Could not determine directory for ' +
                    'smart-r plugin')
        }

        constants.pluginScriptDirectory = new File(smartRDir.canonicalPath, 'HeimScripts')
        log.info("Directory for heim scripts is ${constants.pluginScriptDirectory}")

        if (!skipRScriptsTransfer(config)) {
            def remoteScriptDirectory =  config.smartR.remoteScriptDirectory
            if (!remoteScriptDirectory) {
                remoteScriptDirectory = DEFAULT_REMOTE_RSCRIPTS_DIRECTORY
            }
            remoteScriptDirectory += new Date().getTime()
            constants.remoteScriptDirectoryDir = remoteScriptDirectory
            log.info("Location for R scripts in the Rserve server is ${constants.remoteScriptDirectoryDir}")

            ctx.getBean(RScriptsSynchronizer).start()
        } else {
            log.info('Skipping copying of R script in development mode with local Rserve')
            constants.remoteScriptDirectoryDir = constants.pluginScriptDirectory.absoluteFile
            ctx.getBean(RScriptsSynchronizer).skip()
        }

        if (ctx.containsBean(TRANSMART_EXTENSIONS_REGISTRY_BEAN_NAME)) {
            ctx.getBean(TRANSMART_EXTENSIONS_REGISTRY_BEAN_NAME)
                    .registerAnalysisTabExtension('smartR', '/SmartR/loadScripts', 'addSmartRPanel')
        }

    }

    private boolean skipRScriptsTransfer(config) {
        (!config.RModules.host ||
                config.RModules.host in ['127.0.0.1', '::1', 'localhost']) &&
                Environment.currentEnvironment == Environment.DEVELOPMENT &&
                !config.smartR.alwaysCopyScripts
    }
}
