<script type="text/ng-template" id="variantmap">
<div ng-controller="VariantMapController">


<div id="variantmap-form">
    <div id="sr-variantmap-variantdb-form">
        <div>
            <button ng-click="fetchVariantDB()" ng-disabled="variantDB.invalid || variantDB.running">Fetch VariantDB Data</button>
        </div>

        <div>
            <label for="variantdb-region" >Regions: </label>
            <input id="variantdb-region" type="text" placeholder="13:10000-12500,15:100-5000" ng-model="variantDB.regions"/>
        </div>

        <div>
            <label for="variantdb-genes">Genes: </label>
            <input id="variantdb-genes" type="text" placeholder="ALB,CCR5,CD4,IL10" ng-model="variantDB.selectedGenes"/>
        </div>

        <div>
            <label for="variantdb-serer">Server: </label>
            <input id="variantdb-server" type="text" ng-model="variantDB.server"/>
        </div>

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

        <div id="sr-variantmap-maf-form">
            <p>Frequency:</p>
            <div>
                <label for="sr-variantmap-global-maf">MAF (lower than)</label>
                <input id="sr-variantmap-global-maf" type="number" min="0" max="1" step="0.05" ng-model="variantDB.misc.globalMAF" style="width: 40px;"/>
            </div>
            <div>
                <label for="sr-variantmap-cohort-maf">Group AF (greater than)</label>
                <input id="sr-variantmap-cohort-maf" type="number" min="0" max="1" step="0.05" ng-model="variantDB.misc.cohortMAF" style="width: 40px;"/>
            </div>
        </div>
    </div>

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
            <label for="pdmap-server">Server: </label>
            <input id="pdmap-server" type="text" ng-model="pdmap.server"/>
        </div>
    </div>
</div>
<br/>
<div id="sr-variantmap-messages">
    <p id="error-msgs">{{messages.error}}</p>
    <p id="loading-msgs">Completed Requests: {{messages.finishedRequests}} / {{messages.totalRequests}}</p>
</div>
<br/>
<variantmap-demo data="variantDB.data" show="variantDB.showViz" style="float: left;"></variantmap-demo>
</script>
