<script type="text/ng-template" id="waterfall">
<div ng-controller="WaterfallController as waterfallconn">

    <tab-container>
        %{--========================================================================================================--}%
        %{-- Fetch Data --}%
        %{--========================================================================================================--}%
        <workflow-tab tab-name="Fetch Data" disabled="fetch.disabled">
            <concept-box style="display: inline-block"
                         concept-group="fetch.conceptBoxes.numData"
                         type="LD-numerical"
                         min="1"
                         max="1"
                         label="Numerical Variables"
                         tooltip="Select a continuous variable from the Data Set Explorer Tree and drag it into the box.">
            </concept-box>
            <br/>
            <br/>
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
                <h2>Display Setting:</h2>
                <fieldset class="waterfall-params">
                    <label for="runAnalysis.params.selLowRange">Low Range </label>
                    <select id="runAnalysis.params.selLowRange">
                        <option value="&lt;" selected="selected">&lt;</option>
                        <option value="&lt;=">&lt;=</option>
                    </select>
                    <input id="runAnalysis.params.txtLowRange">

                    <label for="runAnalysis.params.selHighRange">High Range </label>
                    <select id="runAnalysis.params.selHighRange">
                        <option value="&gt;" selected="selected">&gt;</option>
                        <option value="&gt;=">&gt;=</option>
                    </select>
                    <input id="runAnalysis.params.txtHighRange">
                </fieldset>
            </div>
            <hr class="sr-divider">
            <run-button button-name="Create Plot"
                        store-results-in="runAnalysis.scriptResults"
                        script-to-run="run"
                        arguments-to-use="runAnalysis.params"
                        running="runAnalysis.running">
            </run-button>
        </workflow-tab>

    </tab-container>
</div>
</script>
