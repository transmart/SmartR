package heim.rserve

import heim.SmartRRuntimeConstants
import org.rosuda.REngine.Rserve.RConnection

/**
 * Created by piotrzakrzewski on 08/10/15.
 */
class RScriptOutputManager {

    private final RConnection conn
    private final File baseDir
    private final File outputPath

    def RScriptOutputManager(RConnection conn,
                             UUID sessionId,
                             UUID taskId,
                             SmartRRuntimeConstants constants) {
        this.conn = conn
        this.baseDir = constants.baseDir
        this.outputPath = generateOutputPath(sessionId, taskId)
    }

    List<File> downloadFiles() {
        if (!outputPath.mkdirs()) {
            throw new IOException("Failed creating $outputPath")
        }

        List<String> files = listFiles()
        files.collect {
            downloadFromRserve(it)
        }
    }

    private File generateOutputPath(UUID session, UUID taskId) {

        [session, taskId].inject(baseDir)  { File accum, b ->
            new File(accum, b.toString())
        }
    }

    private List<String> listFiles() {
        RUtil.runRCommand(conn, 'list.files()').asNativeJavaObject() as List
    }

    private File downloadFromRserve(String name) {
        File targetFile = new File(outputPath, name)
        if (targetFile.parentFile != outputPath) {
            throw new RuntimeException("$targetFile is a child of $outputPath")
        }

        conn.openFile(name).withStream { InputStream is ->
            targetFile.withOutputStream { OutputStream os ->
                os << is
            }
        }

        targetFile
    }
}
