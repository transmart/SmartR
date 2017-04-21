<script type="text/ng-template" id="rnaseq">
<div ng-controller="RNASeqController as rnaseqconn">

    <tab-container>
        %{--========================================================================================================--}%
        %{-- Fetch Data --}%
        %{--========================================================================================================--}%
        <workflow-tab tab-name="Fetch Data" disabled="fetch.disabled">
            <concept-box style="display: inline-block"
                         concept-group="fetch.conceptBoxes.highDimensional"
                         type="HD"
                         min="1"
                         max="-1"
                         label="RNASeq Data"
                         tooltip="Drag and drop RNASeq data here">
            </concept-box>
            <br/>
            <br/>
            <label>Drag and drop clinical variables to define multiple groups. Please keep in mind that only one type of variables can be compared, e.g. gender (female) with gender (male); not gender (female) with age (>60)!</label>
            <br/>
            <label>Please group the RNASeq data by either categorical variables</label>
            <concept-box style="display: inline-block;"
                         concept-group="fetch.conceptBoxes.cat_grouping"
                         type="LD-categorical"
                         min="-1"
                         max="2"
                         label="Categorical Variables"
                         tooltip="Select two categorical variables to group your data.">
            </concept-box>
            <br/>
            <label>Or by a numerical variable with two different ranges</label>
            <concept-box style="display: inline-block;"
                         concept-group="fetch.conceptBoxes.num_grouping"
                         type="LD-numerical"
                         min="-1"
                         max="2"
                         label="Numerical Variable"
                         tooltip="Select the same numerical variable with two different ranges that you would like to group your data.">
            </concept-box>

            <fetch-button concept-map="fetch.conceptBoxes"
                          loaded="fetch.loaded"
                          running="fetch.running"
                          allowed-cohorts="[1]">
            </fetch-button>
        </workflow-tab>

        %{--========================================================================================================--}%
        %{--Run Analysis--}%
        %{--========================================================================================================--}%
        <workflow-tab tab-name="Run Analysis" disabled="runAnalysis.disabled">
            <div class="heim-input-field sr-input-area">
                <h2>Analysis type:</h2>
                <fieldset class="heim-radiogroup">
                    <label>
                        <input type="radio"
                               ng-model="runAnalysis.params.analysis_type"
                               value="unpaired_group" checked> Two group unpaired
                    </label>
                    <label>
                        <input type="radio"
                               ng-model="runAnalysis.params.analysis_type"
                               value="multi_group"> Multi-group
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
            <rnaseq-plot data="runAnalysis.scriptResults" width="1500" height="1500"></rnaseq-plot>
        </workflow-tab>

    </tab-container>
</div>
</script>
