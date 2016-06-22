package heim.tasks

import au.com.bytecode.opencsv.CSVWriter
import com.google.common.base.Charsets
import com.google.common.base.Stopwatch
import com.google.common.collect.ImmutableMap
import com.google.common.collect.Iterators
import com.google.common.collect.PeekingIterator
import com.google.common.io.Closer
import grails.converters.JSON
import groovy.transform.ToString
import groovy.util.logging.Log4j
import heim.rserve.RServeSession
import heim.rserve.RUtil
import org.codehaus.groovy.grails.support.PersistenceContextInterceptor
import org.rosuda.REngine.REXP
import org.rosuda.REngine.Rserve.RConnection
import org.springframework.beans.factory.annotation.Autowired
import org.springframework.context.annotation.Scope
import org.springframework.stereotype.Component
import org.transmartproject.core.dataquery.DataColumn
import org.transmartproject.core.dataquery.DataRow
import org.transmartproject.core.dataquery.TabularResult
import org.transmartproject.core.dataquery.assay.Assay
import org.transmartproject.core.dataquery.clinical.ClinicalDataResource
import org.transmartproject.core.dataquery.clinical.ClinicalVariable
import org.transmartproject.core.dataquery.clinical.ClinicalVariableColumn
import org.transmartproject.core.dataquery.clinical.PatientRow
import org.transmartproject.core.dataquery.highdim.AssayColumn
import org.transmartproject.core.dataquery.highdim.BioMarkerDataRow
import org.transmartproject.core.dataquery.highdim.HighDimensionDataTypeResource
import org.transmartproject.core.dataquery.highdim.HighDimensionResource
import org.transmartproject.core.dataquery.highdim.assayconstraints.AssayConstraint
import org.transmartproject.core.exceptions.InvalidArgumentsException
import org.transmartproject.core.exceptions.UnexpectedResultException
import org.transmartproject.core.ontology.OntologyTerm
import org.transmartproject.core.querytool.QueriesResource

import java.util.concurrent.*
import java.util.concurrent.atomic.AtomicLong

import static org.transmartproject.core.dataquery.highdim.projections.Projection.DEFAULT_REAL_PROJECTION
import static org.transmartproject.core.ontology.OntologyTerm.VisualAttributes.HIGH_DIMENSIONAL

@Component
@Log4j
@Scope('prototype')
@ToString(includes = 'ontologyTerms,resultInstanceIds,dataType')
class DataFetchTask extends AbstractTask {

    private static final char SEPARATOR = "\t" as char

    private static final int MAX_SIMULTANEOUS_QUERYING = 4
    private static final int MAX_MINUTES_TO_WAIT_FOR_TABRES_FETCH = 2

    private static final AtomicLong THREAD_ID = new AtomicLong()

    Map<String /* label prefix */, OntologyTerm> ontologyTerms

    List<Long> resultInstanceIds /* can have nulls */

    Map<String, Map> assayConstraints

    Map<String, Map> dataConstraints

    String projection

    @Autowired
    private HighDimensionResource highDimensionResource

    @Autowired
    private ClinicalDataResource clinicalDataResource

    @Autowired
    private QueriesResource queriesResource

    @Autowired
    private RServeSession rServeSession // session scoped

    @Autowired
    private PersistenceContextInterceptor interceptor

    private final BlockingQueue queryResultsQueue = new LinkedBlockingQueue()

    private ExecutorService queriesExecutor

    private class GetTabularResultRunnable implements Runnable {

        String label

        Long resultInstanceId

        OntologyTerm ontologyTerm

        @Override
        void run() {
            try {
                if (log.debugEnabled) {
                    log.debug("Will now fetch data with label '$label' for " +
                            "ontology term $ontologyTerm, result instance id " +
                            "$resultInstanceId")
                }

                TabularResult<?, ?> tabularResult = runWithInterceptor this.&doRun

                log.debug("Finished opening tabular result for data set $label")
                queryResultsQueue.put([(label): tabularResult])
            } catch (Exception e) {
                log.error("Fetching of data set $label failed: ${e.message}", e)
                queryResultsQueue.put(e)
            }
        }

        private TabularResult<?, ?> doRun() {
            final Stopwatch stopwatch = Stopwatch.createStarted()

            def res
            if (HIGH_DIMENSIONAL in ontologyTerm.visualAttributes) {
                res = createHighDimensionalResult ontologyTerm, resultInstanceId
            } else {
                res = createClinicalDataResult ontologyTerm, resultInstanceId
            }

            log.info("Fetch for $ontologyTerm (rid $resultInstanceId) " +
                    "finished in $stopwatch")
            res
        }

        private <T> T runWithInterceptor(Closure<T> closure) throws Exception {
            interceptor.init()
            try {
                T ret = closure.call()
                interceptor.flush()
                ret
            } finally {
                interceptor.destroy()
            }
        }
    }

    TaskResult call() throws Exception {
        assert resultInstanceIds.size() > 0 && resultInstanceIds.size() <= 2

        /* actualResultInstanceIds removes nulls, but leaves one element with
         * null if the array would be empty otherwise (note that only in dev
         * is it allowed for resultInstanceIds to be == [null]). */
        def actualResultInstanceIds = resultInstanceIds.findAll() ?: [null]
        List<Map> allDatasets =
                [ontologyTerms.entrySet(), actualResultInstanceIds]
                .combinations()

        queriesExecutor = Executors.newFixedThreadPool(
                Math.min(MAX_SIMULTANEOUS_QUERYING, allDatasets.size()),
                { Runnable r ->
                    new Thread(r).with {
                        daemon = true
                        name = "DataFetchQuery-" + THREAD_ID.incrementAndGet()
                        it
                    }
                } as ThreadFactory)

        clearPreviousCreatedFiles()
        clearPreviousLoadedVariables()

        allDatasets.each {
            Map.Entry<String, OntologyTerm> ontologyTermEntry = it[0]
            Long resultInstanceId = it[1]

            String label = ontologyTermEntry.key + '_s' +
                    (resultInstanceIds.indexOf(resultInstanceId) + 1)
            def runnable = new GetTabularResultRunnable(
                    label: label,
                    resultInstanceId: resultInstanceId,
                    ontologyTerm: ontologyTermEntry.value)

            queriesExecutor.submit(runnable)
        }

        List<String> currentLabels
        int taken = 0
        while (!Thread.currentThread().isInterrupted()
                && taken < allDatasets.size()) {
            def fetchResult = queryResultsQueue.take() // will block
            taken++

            log.debug("Took from queryResultsQueye: $fetchResult")

            if (fetchResult instanceof Exception) {
                log.warn('Error fetching one of the datasets, ' +
                        'aborting writing in the R session (some datasets ' +
                        'may have been written): ' +
                        fetchResult.message, fetchResult)
                return new TaskResult(
                        successful: false,
                        exception: fetchResult)
            }

            assert fetchResult instanceof Map<String, TabularResult>
            try {
                String fileName = writeTabularResult(fetchResult.values().first())
                currentLabels = loadFile(fileName, fetchResult.keySet().first())
            } finally {
                fetchResult.values().first().close()
            }
        }

        if (Thread.interrupted()) {
            throw new InterruptedException("Task was interrupted")
        }

        writeParameters()

        new TaskResult(
                successful: true,
                artifacts: ImmutableMap.of('currentLabels', currentLabels),)
    }

    private clearPreviousCreatedFiles() {
        rServeSession.doWithRConnection { RConnection conn ->
            RUtil.runRCommand conn, 'file.remove(list.files())'
        }
    }

    private clearPreviousLoadedVariables() {
        String removeLoaded = "if (exists('loaded_variables')) { remove(loaded_variables, pos = '.GlobalEnv')}"
        String removePreprocessed = "if (exists('preprocessed')) { remove(preprocessed, pos = '.GlobalEnv')}"
        rServeSession.doWithRConnection { RConnection conn ->
            RUtil.runRCommand conn, removeLoaded
            RUtil.runRCommand conn, removePreprocessed
        }
    }

    private Map<String, Object> packParameters() {
        [
                ontologyTerms: ontologyTerms,
                resultInstanceIds: resultInstanceIds,
                assayConstraints: assayConstraints,
                dataConstraints: dataConstraints,
                projection: projection,
        ]
    }

    private void writeParameters() {
        String loadJsonLite = 'library(jsonlite)'
        def paramsJson = packParameters() as JSON
        String assignParams = "fetch_params <- fromJSON('" +
                RUtil.escapeRStringContent(paramsJson.toString(false)) + "')"

        rServeSession.doWithRConnection { RConnection conn ->
            RUtil.runRCommand conn, loadJsonLite
            RUtil.runRCommand conn, assignParams
        }
    }

    private String /* filename */ writeTabularResult(TabularResult<?, ?> tabularResult) {
        String filename = UUID.randomUUID().toString()
        rServeSession.doWithRConnection { RConnection conn ->
            OutputStream os = new BufferedOutputStream(
                    conn.createFile(filename), 81920)

            log.info("Will start writing tabular result in file $filename")
            final Stopwatch stopwatch = Stopwatch.createStarted()

            Writer writer = new OutputStreamWriter(os, Charsets.UTF_8)

            def csvWriter = new CSVWriter(writer, SEPARATOR)

            try {
                PeekingIterator<? extends DataRow> it =
                        Iterators.peekingIterator(tabularResult.iterator())
                boolean isBioMarker = false 
                try { 
                    isBioMarker = it.peek().hasProperty('bioMarker')
                } catch(NoSuchElementException e) {}
                if (isBioMarker == null) isBioMarker = false
                writeHeader(csvWriter, isBioMarker, tabularResult.indicesList)
                it.each { DataRow row ->
                    if (Thread.interrupted()) {
                        throw new InterruptedException(
                                'Thread was interrupted while dumping the ' +
                                        'TabularResult into a file')
                    }
                    writeLine csvWriter, isBioMarker, row
                }
            } finally {
                csvWriter.close()
            }


            log.info("Finished writing file $filename in $stopwatch")
        }

        filename
    }

    private List<String> /* current labels */ loadFile(String filename, String label) {
        def escapedFilename = RUtil.escapeRStringContent(filename)
        def escapedLabel = RUtil.escapeRStringContent(label)

        List<String> commands = [
                "if (!exists('loaded_variables')) { loaded_variables <- list() }",
                """
                loaded_variables[['$escapedLabel']] <- read.csv(
                               '$escapedFilename', sep = "\t", header = TRUE, stringsAsFactors = FALSE);
                loaded_variables <- loaded_variables[order(names(loaded_variables))]; # for determinism
                names(loaded_variables)""",
        ]
        REXP rexp = rServeSession.doWithRConnection { RConnection conn ->
            RUtil.runRCommand conn, commands[0]
            RUtil.runRCommand conn, commands[1] /* return value */
        }
        rexp.asNativeJavaObject() as List
    }

    private static void writeHeader(CSVWriter writer,
                                    boolean isBioMarker,
                                    List<? extends DataColumn> columns) {
        List line = ['Row Label']
        if (isBioMarker) {
            line += 'Bio marker'
        }
        line += columns.collect {
            if (it instanceof AssayColumn) {
                it.patient.id  // for high dim data, use patient id rather than row label as the CSV header
            } else {
                it.label
            }
        }

        writer.writeNext(line as String[])
    }

    private static void writeLine(CSVWriter writer, boolean isBioMarker, DataRow row) {
        def rowLabel = row.hasProperty('patient') && row.patient.hasProperty('id') ? row.patient.id : row.label
        List line = [rowLabel]
        if (isBioMarker) {
            line += ((BioMarkerDataRow) row).bioMarker
        }
        row.iterator().each {
            line += it
        }

        writer.writeNext(line as String[])
    }

    private TabularResult<AssayColumn, ?> createHighDimensionalResult(
            OntologyTerm ontologyTerm,
            Long resultInstanceId
    ) {
        def subResource = determineSubResource(ontologyTerm)

        def builtAssayConstraints = [
                subResource.createAssayConstraint(
                        AssayConstraint.ONTOLOGY_TERM_CONSTRAINT,
                        concept_key: ontologyTerm.key),
        ]
        if (resultInstanceId) {
            builtAssayConstraints << subResource.createAssayConstraint(
                    AssayConstraint.PATIENT_SET_CONSTRAINT,
                    result_instance_id: resultInstanceId)
        }
        if (assayConstraints) {
            builtAssayConstraints +=
                    assayConstraints.collect { k, v ->
                        subResource.createAssayConstraint(v, k)
                    }
        }

        def builtDataConstraints = []
        if (dataConstraints) {
            builtDataConstraints +=
                    dataConstraints.collect { k, v ->
                        subResource.createDataConstraint(v, k)
                    }
        }

        def builtProjection = projection ?
                subResource.createProjection(projection) :
                subResource.createProjection(DEFAULT_REAL_PROJECTION)

        subResource.retrieveData(
                builtAssayConstraints,
                builtDataConstraints,
                builtProjection)
    }

    private HighDimensionDataTypeResource determineSubResource(OntologyTerm term) {
        Map<HighDimensionDataTypeResource, Collection<Assay>> res =
            highDimensionResource.getSubResourcesAssayMultiMap([
                    highDimensionResource.createAssayConstraint(
                            AssayConstraint.ONTOLOGY_TERM_CONSTRAINT,
                            concept_key: term.key)])

        if (res.keySet().size() > 1) {
            throw new UnexpectedResultException("Found multiple possible " +
                    "data types for ontology term $term: ${res.keySet()}")
        } else if (res.keySet().size() == 0) {
            throw new UnexpectedResultException("No data type found " +
                    "associated with ontology term ${term}.")
        }

        res.keySet().first()
    }

    private TabularResult<ClinicalVariableColumn, PatientRow> createClinicalDataResult(
            OntologyTerm ontologyTerm,
            Long resultInstanceId
    ) {
        if (!resultInstanceId) {
            throw new InvalidArgumentsException(
                    'A result instance id should have been given')
        }

        def clinicalVariable = clinicalDataResource.createClinicalVariable(
                ClinicalVariable.NORMALIZED_LEAFS_VARIABLE,
                concept_path: ontologyTerm.fullName)


        def queryResult = queriesResource.getQueryResultFromId(resultInstanceId)
        clinicalDataResource.retrieveData(
                queryResult,
                [clinicalVariable])
    }

    @Override
    void close() throws Exception {
        // retrieveData() may not work well with interruptions,
        // try to call just shutdown() before shutdownNow()
        try {
            log.debug("Shutting down query executor in $this")
            queriesExecutor.shutdown()
            try {
                def terminated = queriesExecutor.awaitTermination(
                        MAX_MINUTES_TO_WAIT_FOR_TABRES_FETCH, TimeUnit.MINUTES)
                log.debug("Finished waiting for termination in $this. " +
                        "Terminated? $terminated")
                if (!terminated) {
                    log.warn("Tabular result fetching threads did not finish " +
                            "in $MAX_MINUTES_TO_WAIT_FOR_TABRES_FETCH " +
                            "minutes. Will attempt to interrupt them.")
                    queriesExecutor.shutdownNow()
                }
            } catch (InterruptedException ie) {
                log.warn("Interrupted while awaiting for termination of " +
                        "TabularResult fetch in $this")
            }
        } finally {
            Closer closer = Closer.create()
            while (!queryResultsQueue.empty) {
                def o = queryResultsQueue.poll()
                if (o instanceof Map) {
                    TabularResult<?, ?> res = o.values().first()
                    closer.register res
                }
            }
            closer.close()
        }
    }
}
