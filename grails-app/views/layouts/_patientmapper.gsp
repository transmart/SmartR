<script type="text/ng-template" id="patientmapper">

    <div ng-controller="PatientmapperController">

        <tab-container>

            <workflow-tab tab-name="Fetch Data" disabled="fetch.disabled">
                <concept-box style="display: inline-block"
                             concept-group="fetch.conceptBoxes.source"
                             type="LD-categorical"
                             min="1"
                             max="-1"
                             label="Source Variables"
                             tooltip="Select one or more categorical variable which will be used to define the new cohorts.">
                </concept-box>
                <concept-box style="display: inline-block"
                             concept-group="fetch.conceptBoxes.target"
                             type="LD-categorical"
                             min="1"
                             max="-1"
                             label="Target Variables"
                             tooltip="Select one or more categorical variable which will be used to define the new cohorts.">
                </concept-box>
                <br/>
                <br/>
                <fetch-button concept-map="fetch.conceptBoxes"
                              loaded="fetch.loaded"
                              running="fetch.running"
                              allowed-cohorts="[1,2]">
                </fetch-button>
                <run-button button-name="Build Cohort"
                            store-results-in="runAnalysis.scriptResults"
                            script-to-run="run"
                            arguments-to-use="runAnalysis.params"
                            running="runAnalysis.running">
                </run-button>
            </workflow-tab>

        </tab-container>

    </div>

</script>
