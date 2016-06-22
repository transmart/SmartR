
<script type="text/ng-template" id="correlation">

    <div ng-controller="CorrelationController">

        <tab-container>

            <workflow-tab tab-name="Fetch Data" disabled="fetch.disabled">
                <concept-box style="display: inline-block"
                             concept-group="fetch.conceptBoxes.datapoints"
                             type="LD-numerical"
                             min="2"
                             max="2"
                             label="Numerical Variables"
                             tooltip="Select two numerical variables from the tree to compare them.">
                </concept-box>
                <concept-box style="display: inline-block;"
                             concept-group="fetch.conceptBoxes.annotations"
                             type="LD-categorical"
                             min="0"
                             max="-1"
                             label="(optional) Categorical Variables"
                             tooltip="Select an arbitrary amount of categorical variables from the tree to annotate the numerical datapoints.">
                </concept-box>
                <br/>
                <br/>
                <fetch-button concept-map="fetch.conceptBoxes"
                              loaded="fetch.loaded"
                              running="fetch.running"
                              allowed-cohorts="[1]">
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
                                   value="log2" checked> Log2
                        </label>
                        <label>
                            <input type="radio"
                                   ng-model="runAnalysis.params.transformation"
                                   value="log10" checked> Log10

                        </label>
                    </fieldset>
                </div>
                <div class="heim-input-field sr-input-area">
                    <h2>Correlation computation method:</h2>
                    <fieldset class="heim-radiogroup">
                        <label>
                            <input type="radio"
                                   ng-model="runAnalysis.params.method"
                                   value="pearson" checked> Pearson
                        </label>
                        <label>
                            <input type="radio"
                                   ng-model="runAnalysis.params.method"
                                   value="kendall"> Kendall
                        </label>
                        <label>
                            <input type="radio"
                                   ng-model="runAnalysis.params.method"
                                   value="spearman"> Spearman
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
                <correlation-plot data="runAnalysis.scriptResults" width="1500" height="1500"></correlation-plot>
            </workflow-tab>

        </tab-container>

    </div>

</script>
