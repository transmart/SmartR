
<script type="text/ng-template" id="correlation">

    <div ng-controller="CorrelationController">

        <tab-container>

            <workflow-tab tab-name="Fetch Data">
                <concept-box style="display: inline-block"
                             concept-group="conceptBoxes.datapoints"
                             type="LD-numerical"
                             min="2"
                             max="2"
                             label="Numerical Variables"
                             tooltip="Select two numerical variables from the tree to compare them.">
                </concept-box>
                <concept-box style="display: inline-block;"
                             concept-group="conceptBoxes.annotations"
                             type="LD-categorical"
                             min="0"
                             max="-1"
                             label="(optional) Categorical Variables"
                             tooltip="Select an arbitrary amount of categorical variables from the tree to annotate the numerical datapoints.">
                </concept-box>
                <br/>
                <br/>
                <fetch-button concept-map="conceptBoxes"
                              show-summary-stats="false"
                                disabled="false">
                </fetch-button>
            </workflow-tab>

            <workflow-tab tab-name="Run Analysis">
                <select ng-model="params.method">
                    <option value="pearson">Pearson (Default)</option>
                    <option value="kendall">Kendall</option>
                    <option value="spearman">Spearman</option>
                </select>

                <br/>
                <br/>
                <run-button button-name="Create Plot"
                            store-results-in="scriptResults"
                            script-to-run="run"
                            arguments-to-use="params"
                            serialized="false"
                            disabled="false">
                </run-button>
                <br/>
                <br/>
                <correlation-plot data="scriptResults" width="1200" height="1200"></correlation-plot>
            </workflow-tab>

        </tab-container>

    </div>

</script>
