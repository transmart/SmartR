
<script type="text/ng-template" id="ppmidemo">
<div ng-controller="PPMIDemoController">


<div id="ppmi-form">
    <div id="sr-ppmi-variantdb-form">
        <div>
            <button ng-click="fetchVariantDB()" ng-disabled="variantDB.invalid">Fetch VariantDB Data</button>
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

        <div id="sr-ppmi-region-form">
            <p>Gene Location:</p>
            <div>
                <input id="sr-ppmi-exonic-check" type="checkbox" />
                <label for="sr-ppmi-exonic-check">exonic</label>
            </div>

            <div>
                <input id="sr-ppmi-intronic-check" type="checkbox" />
                <label for="sr-ppmi-intronic-check">intronic</label>
            </div>

            <div>
                <input id="sr-ppmi-upstream-check" type="checkbox" />
                <label for="sr-ppmi-upstream-check">upstream</label>
            </div>

            <div>
                <input id="sr-ppmi-downstream-check" type="checkbox" />
                <label for="sr-ppmi-downstream-check">downstream</label>
            </div>

            <div>
                <input id="sr-ppmi-splicing-check" type="checkbox" />
                <label for="sr-ppmi-slicing-check">slicing</label>
            </div>

            <div>
                <input id="sr-ppmi-intergenetic-check" type="checkbox" />
                <label for="sr-ppmi-intergenetic-check">intergenetic</label>
            </div>
        </div>

        <div id="sr-ppmi-exonic-function-form">
            <p>Exonic Function</p>
            <div>
                <input id="sr-ppmi-non-synonymous-check" type="checkbox" />
                <label for="sr-ppmi-non-synonymous-check">non-synonymous SNV</label>
            </div>

            <div>
                <input id="sr-ppmi-synonymous-check" type="checkbox" />
                <label for="sr-ppmi-synonymous-check">synonymous SNV</label>
            </div>

            <div>
                <input id="sr-ppmi-insertion-check" type="checkbox" />
                <label for="sr-ppmi-insertion-check">insertion</label>
            </div>

            <div>
                <input id="sr-ppmi-deletion-check" type="checkbox" />
                <label for="sr-ppmi-deletion-check">deletion</label>
            </div>

            <div>
                <input id="sr-ppmi-block-substitution-check" type="checkbox" />
                <label for="sr-ppmi-block-substitution-check">block-substitution</label>
            </div>
        </div>

        <div id="sr-ppmi-maf-form">
            <p>Misc:</p>
            <label for="sr-ppmi-maf">M.A.F.:</label>
            <input id="sr-ppmi-maf" type="number" min="0" max="1" value="0.1" step="0.05" style="width: 40px;"/>
        </div>
    </div>

    <div id="sr-ppmi-pdmap-form">
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
<div id="sr-ppmi-messages">
    <p id="error-msgs">{{messages.error}}</p>
    <p id="loading-msgs">Completed Requests: {{messages.finishedRequests}} / {{messages.totalRequests}}</p>
</div>
</script>
