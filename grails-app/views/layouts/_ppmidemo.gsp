
<script type="text/ng-template" id="ppmidemo">
<div ng-controller="PPMIDemoController">


<div id="ppmi-form">
    <div>
        <label for="variantdb-region" >Regions: </label>
        <input id="variantdb-region" type="text" placeholder="13:10000-12500,15:100-5000" ng-model="variantDB.regions"/>

        <label for="variantdb-genes">Genes: </label>
        <input id="variantdb-genes" type="text" placeholder="ALB,CCR5,CD4,IL10" ng-model="variantDB.genes"/>

        <label for="variantdb-serer">Server: </label>
        <input id="variantdb-server" type="text" ng-model="variantDB.server"/>
    </div>

    <div>
        <button ng-click="fetchVariantDB()" ng-disabled="variantDB.invalid">Fetch VariantDB Data</button>
    </div>

    <div>
        <label for="pdmap-user">UserID: </label>
        <input id="pdmap-user" type="text" ng-model="pdmap.user"/>

        <label for="pdmap-password">Password: </label>
        <input id="pdmap-password" type="password" ng-model="pdmap.password"/>

        <label for="pdmap-server">Server: </label>
        <input id="pdmap-server" type="text" ng-model="pdmap.server"/>
    </div>

    <div>
        <button ng-click="createPDMapLayout()" ng-disabled="pdmap.invalid">Create PDMap Overlay</button>
    </div>
</div>

<p id="error-msgs">{{messages.error}}</p>

</div>
</script>
