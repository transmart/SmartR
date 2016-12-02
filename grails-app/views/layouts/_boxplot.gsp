
<script type="text/ng-template" id="boxplot">

<div ng-controller="BoxplotController">

    <tab-container>

        <workflow-tab tab-name="Fetch Data" disabled="fetch.disabled">
            <concept-box style="display: inline-block;"
                         concept-group="fetch.conceptBoxes.highDimensional"
                         type="HD"
                         min="-1"
                         max="-1"
                         label="High Dimensional Variables"
                         tooltip="Select one or more high dimensional variables that you would like to have displayed.">
            </concept-box>
            <concept-box style="display: inline-block;"
                         concept-group="fetch.conceptBoxes.numData"
                         type="LD-numerical"
                         min="-1"
                         max="-1"
                         label="Numerical Variables"
                         tooltip="Select one or more numerical variables that you would like to have displayed.">
            </concept-box>
            <concept-box style="display: inline-block;"
                         concept-group="fetch.conceptBoxes.groups"
                         type="LD-categorical"
                         min="-1"
                         max="-1"
                         label="(optional) Categorical Variables"
                         tooltip="Select one or more categorical variables to group your data.">
            </concept-box>
            <br/>
            <br/>
            <biomarker-selection biomarkers="fetch.selectedBiomarkers"></biomarker-selection>
            <hr class="sr-divider">
            <fetch-button concept-map="fetch.conceptBoxes"
                          loaded="fetch.loaded"
                          running="fetch.running"
                          biomarkers="fetch.selectedBiomarkers"
                          disabled="fetch.button.disabled"
                          message="fetch.button.message"
                          allowed-cohorts="[1,2]">
            </fetch-button>
        </workflow-tab>

        <workflow-tab tab-name="Run Analysis" disabled="runAnalysis.disabled">
           <div class="heim-input-field sr-input-area">
                <h2>Data transformation:</h2>
                <fieldset class="heim-radiogroup">
                    <label>
                        <input type="radio"
                               ng-model="runAnalysis.params.transformation"
                               value="raw" checked> Raw Values
                    </label>
                    <label>
                        <input type="radio"
                               ng-model="runAnalysis.params.transformation"
                               value="log2"> Log2
                    </label>
                    <label>
                        <input type="radio"
                               ng-model="runAnalysis.params.transformation"
                               value="log10"> Log10

                    </label>
                </fieldset>
            </div>
            <hr class="sr-divider">
            <run-button button-name="Create Plot"
                        store-results-in="runAnalysis.scriptResults"
                        script-to-run="run"
                        arguments-to-use="runAnalysis.params"
                        running="runAnalysis.running">
            </run-button>
            <br/>
            <br/>
            <boxplot data="runAnalysis.scriptResults"></boxplot>
        </workflow-tab>

    </tab-container>

</div>

</script>
