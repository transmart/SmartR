
<script type="text/ng-template" id="boxplot">

<div ng-controller="BoxplotController">

    <tab-container>

        <workflow-tab tab-name="Fetch Data" disabled="fetch.disabled">
            <concept-box style="display: inline-block;"
                         concept-group="fetch.conceptBoxes.datapoints"
                         type="LD-numerical"
                         min="1"
                         max="1"
                         label="Numerical Variable"
                         tooltip="Select a single numerical variable that you would like to have displayed.">
            </concept-box>
            %{--Nice idea but somehow lost it's initial purpose because cross-study support is gone.
                Maybe implement later--}%
            %{--<concept-box style="display: inline-block;"--}%
                         %{--concept-group="fetch.conceptBoxes.subsets"--}%
                         %{--type="LD-categorical"--}%
                         %{--min="0"--}%
                         %{--max="-1"--}%
                         %{--label="(optional) Categorical Variables"--}%
                         %{--tooltip="Select an arbitrary amount of categorical variables to split the numerical variable into subsets.">--}%
            %{--</concept-box>--}%
            <br/>
            <br/>
            <fetch-button concept-map="fetch.conceptBoxes"
                          loaded="fetch.loaded"
                          running="fetch.running">
            </fetch-button>
        </workflow-tab>

        <workflow-tab tab-name="Run Analysis" disabled="runAnalysis.disabled">
            <br/>
            <br/>
            <run-button button-name="Create Plot"
                        store-results-in="runAnalysis.scriptResults"
                        script-to-run="run"
                        arguments-to-use="runAnalysis.params"
                        running="runAnalysis.running">
            </run-button>
            <br/>
            <br/>
            <boxplot data="runAnalysis.scriptResults" width="1000" height="500"></boxplot>
        </workflow-tab>

    </tab-container>

</div>

</script>
