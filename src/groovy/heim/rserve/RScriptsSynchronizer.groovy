package heim.rserve

import com.google.common.base.Charsets
import groovy.transform.TypeChecked
import groovy.util.logging.Log4j
import heim.SmartRRuntimeConstants
import org.rosuda.REngine.Rserve.RConnection
import org.rosuda.REngine.Rserve.RFileOutputStream
import org.springframework.beans.factory.annotation.Autowired
import org.springframework.stereotype.Component

import java.nio.file.FileVisitResult
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.SimpleFileVisitor
import java.nio.file.attribute.BasicFileAttributes
import java.util.zip.ZipEntry
import java.util.zip.ZipOutputStream

/**
 * Synchronizes the R script from the machine running the web application
 * to the machine running Rserve.
 */
@Component
@TypeChecked
@Log4j
class RScriptsSynchronizer {

    @Autowired
    SmartRRuntimeConstants constants

    public static final String BUNDLE_FIX = 'smartR_bundle.zip'

    private static final int MAX_NUMBER_OF_RETRIES = 10
    private static final int RETRY_INTERVAL = 60L // seconds

    @Autowired
    private RConnectionProvider rConnectionProvider

    private boolean copySuccessful = false
    private final Object lock = new Object()

    void start() {
        Thread.start('r-scripts-synchronizer') {
            int totalTries = MAX_NUMBER_OF_RETRIES + 1
            while (totalTries-- > 0) {
                try {
                    ensureScriptsAreCopied()
                    return
                } catch (Throwable t) {
                    if (t instanceof InterruptedException) {
                        log.warn('R scripts synchronization interrupted')
                        return
                    }
                    log.warn('Synchronization of R scripts has failed. ' +
                            'Retries left: ' + totalTries, t)
                    Thread.sleep(RETRY_INTERVAL * 1000L)
                }
            }

            log.error('Synchronization of R scripts has failed. No more retries')
        }
    }

    void ensureScriptsAreCopied() {
        synchronized (lock) {
            if(!copySuccessful) {
                try {
                    copyFilesOver()
                } catch (Throwable t) {
                    throw new IllegalStateException("Can't copy the smartR scripts to the Rserver machine", t)
                }
                copySuccessful = true
                log.info('Finished copying over the smartR scripts to the Rserve machine')
            }
        }
    }

    void skip() {
        synchronized (lock) {
            copySuccessful = true
            log.info('Skipping on copying over the smartR scripts to the Rserve machine')
        }
    }

    private void copyFilesOver() {
        RConnection rConnection = rConnectionProvider.get()
        try {
            RFileOutputStream os = rConnection.createFile(BUNDLE_FIX)
            // buffering the output improves performance (o/wise writing just
            // one byte would send a packet on itself) and goes around a bug
            // in RFileOutputStream::write(byte b), which always writes a NUL
            // byte, no matter the argument.
            BufferedOutputStream bos = new BufferedOutputStream(os)
            createZip(bos)

            String escapedRemoteDir =
                    RUtil.escapeRStringContent(constants.remoteScriptDirectoryDir.absoluteFile.toString())
            assert escapedRemoteDir != '/'
            RUtil.runRCommand(rConnection,
                    "unlink('$escapedRemoteDir', recursive = TRUE)")

            String escapedZipName = RUtil.escapeRStringContent(BUNDLE_FIX)
            RUtil.runRCommand(rConnection,
                    "unzip('$escapedZipName', exdir = '$escapedRemoteDir')")
        } finally {
            rConnection.close()
        }
    }

    /* the passed in OutputStrem will be closed */
    private void createZip(OutputStream destinationStream) {
        ZipOutputStream zs = new ZipOutputStream(destinationStream, Charsets.UTF_8)
        try {
            Path basePath = constants.pluginScriptDirectory.toPath()

            Files.walkFileTree(
                    basePath,
                    new RScriptVisitor(zs: zs, basePath: basePath))
        } finally {
            if (zs) {
                zs.close()
            } else {
                destinationStream.close()
            }
        }
    }
    private static class RScriptVisitor extends SimpleFileVisitor<Path> {
        ZipOutputStream zs
        Path basePath

        @Override
        FileVisitResult visitFile(Path file, BasicFileAttributes attrs) throws IOException {
            Path relativePath = basePath.relativize(file)
            ZipEntry entry = new ZipEntry(relativePath.toString())

            zs.putNextEntry(entry)
            Files.copy(file, zs)
            zs.closeEntry()

            super.visitFile(file, attrs)
        }
    }
}
