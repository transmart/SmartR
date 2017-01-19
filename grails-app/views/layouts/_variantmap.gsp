<script type="text/ng-template" id="variantmap">
<div ng-controller="VariantMapController">

<tab-container>
    <workflow-tab tab-name="Fetch Data" disabled="fetch.disabled">
        <concept-box
            concept-group="fetch.conceptBoxes.highDimensional"
            type="HD"
            min="0"
            max="-1"
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
            </br>
            <div id="sr-variantmap-region-form">
                <p>Gene Location:</p>
                <div>
                    <input id="sr-variantmap-exonic-check" type="checkbox" ng-model="variantDB.func_refgene.exonic"/>
                    <label for="sr-variantmap-exonic-check">exonic</label>
                </div>

                <div>
                    <input id="sr-variantmap-intronic-check" type="checkbox" ng-model="variantDB.func_refgene.intronic"/>
                    <label for="sr-variantmap-intronic-check">intronic</label>
                </div>

                <div>
                    <input id="sr-variantmap-upstream-check" type="checkbox" ng-model="variantDB.func_refgene.upstream"/>
                    <label for="sr-variantmap-upstream-check">upstream</label>
                </div>

                <div>
                    <input id="sr-variantmap-downstream-check" type="checkbox" ng-model="variantDB.func_refgene.downstream"/>
                    <label for="sr-variantmap-downstream-check">downstream</label>
                </div>

                <div>
                    <input id="sr-variantmap-splicing-check" type="checkbox" ng-model="variantDB.func_refgene.splicing"/>
                    <label for="sr-variantmap-splicing-check">splicing</label>
                </div>

                <div>
                    <input id="sr-variantmap-intergenetic-check" type="checkbox" ng-model="variantDB.func_refgene.intergenetic"/>
                    <label for="sr-variantmap-intergenetic-check">intergenetic</label>
                </div>
            </div>
            </br>
            <div id="sr-variantmap-exonic-function-form">
                <p>Exonic Function</p>
                <div>
                    <input id="sr-variantmap-synonymous-check" type="checkbox" ng-model="variantDB.exonicfunc_refgene.synonymous_SNV"/>
                    <label for="sr-variantmap-synonymous-check">synonymous SNV</label>
                </div>

                <div>
                    <input id="sr-variantmap-non-synonymous-check" type="checkbox" ng-model="variantDB.exonicfunc_refgene.nonsynonymous_SNV"/>
                    <label for="sr-variantmap-non-synonymous-check">non-synonymous SNV</label>
                </div>

                <div>
                    <input id="sr-variantmap-insertion-check" type="checkbox" ng-model="variantDB.exonicfunc_refgene.frameshift_insertion"/>
                    <label for="sr-variantmap-insertion-check">frameshift insertion</label>
                </div>

                <div>
                    <input id="sr-variantmap-insertion-check" type="checkbox" ng-model="variantDB.exonicfunc_refgene.nonframeshift_insertion"/>
                    <label for="sr-variantmap-insertion-check">nonframeshift insertion</label>
                </div>

                <div>
                    <input id="sr-variantmap-deletion-check" type="checkbox" ng-model="variantDB.exonicfunc_refgene.frameshift_deletion"/>
                    <label for="sr-variantmap-deletion-check">frameshift deletion</label>
                </div>

                <div>
                    <input id="sr-variantmap-deletion-check" type="checkbox" ng-model="variantDB.exonicfunc_refgene.nonframeshift_deletion"/>
                    <label for="sr-variantmap-deletion-check">nonframeshift deletion</label>
                </div>

                <div>
                    <input id="sr-variantmap-block-substitution-check" type="checkbox" ng-model="variantDB.exonicfunc_refgene.frameshift_substitution"/>
                    <label for="sr-variantmap-block-substitution-check">frameshift substitution</label>
                </div>

                <div>
                    <input id="sr-variantmap-block-substitution-check" type="checkbox" ng-model="variantDB.exonicfunc_refgene.nonframeshift_substitution"/>
                    <label for="sr-variantmap-block-substitution-check">nonframeshift substitution</label>
                </div>
            </div>
            </br>
            <div id="sr-variantmap-maf-form">
                <p>Frequency:</p>
                <div>
                    <label for="sr-variantmap-global-maf">MAF (lower than)</label>
                    <input id="sr-variantmap-global-maf" type="number" min="0.2" max="1" step="0.05" ng-model="variantDB.misc.globalMAF" style="width: 40px;"/>
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
        <div id="sr-variantmap-pdmap-form">
            <div>
                <button ng-click="createPDMapLayout()" ng-disabled="pdmap.invalid">Create PDMap Overlay</button>
            </div>

            <div>
                <label for="pdmap-user">UserID: </label>
                <input id="pdmap-user" type="text" ng-model="pdmap.user"/>
            </div>

            <div>
                <label for="pdmap-password">Password: </label>
                <input id="pdmap-password" type="password" ng-model="pdmap.password"/>
            </div>

            <div>
                <label for="pdmap-server">VariantDB Endpoint: </label>
                <input id="pdmap-server" type="text" ng-model="pdmap.server"/>
            </div>
        </div>
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
