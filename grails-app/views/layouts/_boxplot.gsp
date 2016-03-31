
<script type="text/ng-template" id="boxplot">

<div ng-controller="BoxplotController">

    <tab-container>

        <workflow-tab tab-name="Fetch Data">
            <concept-box style="display: inline-block;"
                         concept-group="conceptBoxes.datapoints"
                         type="LD-numerical"
                         min="1"
                         max="1"
                         label="Numerical Variable"
                         tooltip="Select a single numerical variable that you would like to have displayed.">
            </concept-box>
            <concept-box style="display: inline-block;"
                         concept-group="conceptBoxes.subsets"
                         type="LD-categorical"
                         min="0"
                         max="-1"
                         label="(optional) Categorical Variables"
                         tooltip="Select an arbitrary amount of categorical variables to split the numerical variable into subsets.">
            </concept-box>
            <br/>
            <br/>
            <fetch-button concept-map="conceptBoxes" show-summary-stats="false" disabled="false"></fetch-button>
        </workflow-tab>

        <workflow-tab tab-name="Run Analysis">
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
            <boxplot data="scriptResults" width="1000" height="500"></boxplot>
        </workflow-tab>

    </tab-container>

</div>

</script>
