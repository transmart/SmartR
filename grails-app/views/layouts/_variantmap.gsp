<script type="text/ng-template" id="variantmap">
<div ng-controller="VariantMapController">

<tab-container>
    <workflow-tab tab-name="Fetch Data" disabled="fetch.disabled">
        <concept-box
            concept-group="fetch.conceptBoxes.highDimensional"
            type="HD"
            min="0"
            max="1"
            label="(optional) High Dimensional Variables"
            tooltip="Select high dimensional data node(s) from the data tree and drag it into the box.
            The nodes need to be from the same platform.">
        </concept-box>

        <concept-box
            concept-group="fetch.conceptBoxes.numeric"
            type="LD-numerical"
            min="0"
            max="-1"
            label="(optional) Numerical Variables"
            tooltip="Select numeric data node(s) from the data tree and drag it into the box.">
        </concept-box>

        <concept-box
            concept-group="fetch.conceptBoxes.categoric"
            type="LD-categorical"
            min="0"
            max="-1"
            label="(optional) Categoric Variables"
            tooltip="Select categoric data node(s) from the data tree and drag it into the box.">
        </concept-box>
        <biomarker-selection biomarkers="fetch.selectedBiomarkers"></biomarker-selection>
        <div id="variantmap-form">
            <div id="sr-variantmap-variantdb-form">
                <div>
                    <label for="variantdb-serer">Server: </label>
                    <input id="variantdb-server" type="text" ng-model="variantDB.server"/>
                </div>
            </div>
        </div>

        <hr class="sr-divider">
        <fetch-button
                callback="fetchVariantDB"
                loaded="fetch.loaded"
                running="fetch.running"
                concept-map="fetch.conceptBoxes"
                biomarkers="fetch.selectedBiomarkers"
                allowed-cohorts="[1,2]"
                has-preprocess-tab="false"
                ng-disabled="!fetch.selectedBiomarkers.length">
        </fetch-button>
        </br>
        <div id="sr-variantmap-messages">
            <p id="error-msgs">{{messages.error}}</p>
        </div>
        <br/>
    </workflow-tab>
    <workflow-tab tab-name="Run Analysis" disabled="runAnalysis.disabled">
        <run-button button-name="Create Plot"
                    store-results-in="runAnalysis.scriptResults"
                    script-to-run="run"
                    arguments-to-use="runAnalysis.params"
                    running="runAnalysis.running">
        </run-button>
        <capture-plot-button filename="variantmap.svg" target="variant-map"></capture-plot-button>
        <br/>
        <br/>
        <variant-map data="runAnalysis.scriptResults"></variant-map>
    </workflow-tab>
</tab-container>
</div>
</script>
