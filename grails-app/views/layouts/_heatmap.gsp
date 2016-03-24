<script type="text/ng-template" id="heatmap">
<div ng-controller="HeatmapController">

    <tab-container>
        %{--========================================================================================================--}%
        %{-- Fetch Data --}%
        %{--========================================================================================================--}%
        <workflow-tab tab-name="Fetch Data">

            <table>
                <tr>
                    <td>
                        <concept-box concept-group="fetch.conceptBoxes.highDimensional"
                                     type="HD"
                                     min="1"
                                     max="-1"
                                     label="High Dimensional Variables"
                                     tooltip="Select high dimensional data node(s) from the Data Set Explorer Tree and drag it into the box. The nodes needs to be from the same platform.">
                        </concept-box>
                    </td>
                    %{--<td>--}%
                        %{--<concept-box concept-group="fetch.conceptBoxes.numerical"--}%
                                     %{--type="LD-numerical"--}%
                                     %{--min="0"--}%
                                     %{--max="-1"--}%
                                     %{--label="(optional) Numerical Variables"--}%
                                     %{--tooltip="Select an arbitrary number of numerical variables to expand the heatmap with low dimensional data.">--}%
                        %{--</concept-box>--}%
                    %{--</td>--}%
                    %{--<td>--}%
                        %{--<concept-box concept-group="fetch.conceptBoxes.categorical"--}%
                                     %{--type="LD-categorical"--}%
                                     %{--min="0"--}%
                                     %{--max="-1"--}%
                                     %{--label="(optional) Categorical Variables"--}%
                                     %{--tooltip="Select an arbitrary number of categorical variables to expand the heatmap with low dimensional data.">--}%
                        %{--</concept-box>--}%
                    %{--</td>--}%
                </tr>
            </table>


            <biomarker-selection biomarkers="fetch.selectedBiomarkers"></biomarker-selection>
            <hr class="sr-divider">
            <fetch-button
                    disabled="fetch.disabled"
                    concept-map="fetch.conceptBoxes"
                    biomarkers="fetch.selectedBiomarkers"
                    show-summary-stats="true"
                    summary-data="fetch.summaryData">
            </fetch-button>
            <br/>
            <summary-stats summary-data="fetch.summaryData"></summary-stats>
        </workflow-tab>

        %{--========================================================================================================--}%
        %{-- Preprocess Data --}%
        %{--========================================================================================================--}%
        <workflow-tab tab-name="Preprocess">
            %{--Aggregate Probes--}%
            <div class="heim-input-field">
                <input type="checkbox" ng-model="preprocess.params.aggregateProbes">
                <span>Aggregate probes</span>
            </div>

            <hr class="sr-divider">

            <preprocess-button params="preprocess.params"
                               show-summary-stats="true"
                               summary-data="preprocess.summaryData"
                               disabled="preprocess.disabled">
            </preprocess-button>

            <br/>
            <summary-stats summary-data="preprocess.summaryData"></summary-stats>
        </workflow-tab>


        %{--========================================================================================================--}%
        %{--Run Analysis--}%
        %{--========================================================================================================--}%
        <workflow-tab tab-name="Run Analysis">
            %{--Number of max row to display--}%
            <div class="heim-input-field heim-input-number sr-input-area">
                Show <input type="text" id="txtMaxRow" ng-model="runAnalysis.params.max_row">
                of {{fetch.summaryData.summary[0].$$state.value.json.data[0].totalNumberOfValuesIncludingMissing /
                    fetch.summaryData.summary[0].$$state.value.json.data[0].numberOfSamples}}
                rows in total. (< 1000 is preferable.)
            </div>

            %{--Type of sorting to apply--}%
            <div class="heim-input-field sr-input-area">
                <h2>Sort on:</h2>
                <fieldset class="heim-radiogroup">
                    <label>
                        <input type="radio" ng-model="runAnalysis.params.sorting" name="sortingSelect" value="nodes"
                               checked> Nodes
                    </label>
                    <label>
                        <input type="radio" ng-model="runAnalysis.params.sorting" name="sortingSelect" value="subjects">
                        Subjects
                    </label>
                </fieldset>
            </div>

            %{--Link to Biocompendium or PDMap--}%
            <div class="heim-input-field  sr-input-area">
                <pd-map-login criteria="runAnalysis.params.ranking"
                            subsets="runAnalysis.subsets"
                            login="runAnalysis.params.pdMapLogin"
                            password="runAnalysis.params.pdMapPassword"
                            pdMapLinkEnabled="runAnalysis.params.pdMapLinkEnabled">
                </pd-map-login>
            </div>

            %{--Type of sorting to apply--}%
            <div class="heim-input-field  sr-input-area">
                <sorting-criteria criteria="runAnalysis.params.ranking" subsets="runAnalysis.subsets">
                </sorting-criteria>
            </div>

            <hr class="sr-divider">

            <run-button button-name="Create Plot"
                        store-results-in="runAnalysis.scriptResults"
                        script-to-run="run"
                        arguments-to-use="runAnalysis.params"
                        serialized="true"
                        disabled="runAnalysis.disabled">
            </run-button>
            <capture-plot-button filename="heatmap.svg" disabled="runAnalysis.download.disabled"></capture-plot-button>
            <download-results-button disabled="runAnalysis.download.disabled"></download-results-button>
            <br/>
            <heatmap-plot data="runAnalysis.scriptResults" width="1200" height="1200"></heatmap-plot>

        </workflow-tab>

    </tab-container>
</div>
</script>
