
<script type="text/ng-template" id="logit">

    <div ng-controller="LogitController">

        <tab-container>

            <workflow-tab tab-name="Fetch Data" disabled="fetch.disabled">
                <concept-box style="display: inline-block"
                             concept-group="fetch.conceptBoxes.datapoints"
                             type="LD-numerical"
                             min="2"
                             max="2"
                             label="Numerical Variables"
                             tooltip="Select two numerical variables from the tree for plotting.">
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
                    <h2>X Data transformation:</h2>
                    <fieldset class="heim-radiogroup">
                        <label>
                            <input type="radio"
                                   ng-model="runAnalysis.params.transformationx"
                                   value="raw" checked> Raw Values
                        </label>
                        <label>
                            <input type="radio"
                                   ng-model="runAnalysis.params.transformationx"
                                   value="log2" checked> Log2
                        </label>
                        <label>
                            <input type="radio"
                                   ng-model="runAnalysis.params.transformationx"
                                   value="log10" checked> Log10

                        </label>
                    </fieldset>
                </div>
                <div class="heim-input-field sr-input-area">
                    <h2>Y Data Transformation:</h2>
                    <fieldset class="heim-radiogroup">
                        <label>
                            <input type="radio"
                                   ng-model="runAnalysis.params.transformationy"
                                   value="raw" checked> Raw Values
                        </label>
                        <label>
                            <input type="radio"
                                   ng-model="runAnalysis.params.transformationy"
                                   value="normalize" checked> normalize
                        </label>
                        <label>
                            <input type="radio"
                                   ng-model="runAnalysis.params.transformationy"
                                   value="clamp" checked> clamp
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
                <capture-plot-button filename="logit.svg" target="logit-plot"></capture-plot-button>
                <br/>
                <br/>
                <logit-plot data="runAnalysis.scriptResults" width="1500" height="1500"></logit-plot>
            </workflow-tab>

        </tab-container>

    </div>

</script>
