//# sourceURL=heatmap.js

function HeatmapController ($scope, smartRUtils, commonWorkflowService) {

    commonWorkflowService.initializeWorkflow('heatmap', $scope);

    // ------------------------------------------------------------- //
    // Fetch data                                                    //
    // ------------------------------------------------------------- //
    $scope.fetch = {
        btn : {
            disabled : true
        },
        conceptBoxes : {
            highDimensional : [],
            numerical : [],
            categorical : []
        },
        selectedBiomarkers : [],
        scriptResults : {},
        totalSamples : 0,
        subsets : 0
    };

    // ------------------------------------------------------------- //
    // Preprocess                                                    //
    // ------------------------------------------------------------- //
    $scope.preprocess = {
        btn : {
            disabled : true
        },
        params :  {
            aggregateProbes : false
        },
        scriptResults : {}
    };

    // ------------------------------------------------------------- //
    // Run Heatmap                                                   //
    // ------------------------------------------------------------- //
    $scope.runAnalysis = {
        btn : {
            disabled : true
        },
        params: {
            max_row : 100,
            sorting : 'nodes',
            ranking : 'coef'
        },
        download : {
            disabled : true
        },
        subsets : 0,
        scriptResults : {}
    };

    commonWorkflowService.registerCondition(
        ['fetch.conceptBoxes.highDimensional'],
        ['fetch.btn'],
        this.whenSelectHighDimensionalNodes
    );

    commonWorkflowService.registerCondition(
        ['fetch.btn.disabled', 'fetch.totalSamples', 'fetch.subsets'],
        ['preprocess', 'runAnalysis'],
        this.whenFetchHeatmapData
    );

    commonWorkflowService.registerCondition(
        ['preprocess.btn.disabled'],
        ['fetch', 'runAnalysis'],
        this.whenPreprocessHeatmapData
    );

    commonWorkflowService.registerCondition(
        ['runAnalysis.btn.disabled', 'runAnalysis.scriptResults'],
        ['fetch', 'preprocess', 'runAnalysis.download'],
        this.whenRunHeatmapAnalysis
    );
}

/**
 * Components are enabled/disabled according to presence of high dimensional nodes
 * @param newValues
 * @param oldValues
 * @param scope
 * @param affectedComponents
 */
HeatmapController.prototype.whenSelectHighDimensionalNodes = function (newValues, oldValues, scope, affectedComponents) {
    affectedComponents.forEach(function (component) {
        if (component.hasOwnProperty('disabled')) {
            component.disabled = newValues[0].length <= 0;
        }
    });
};

/**
 * Collective behaviour when fetching data for heatmap
 * @param newValues - new values of components models that determine other components states
 * @param oldValues - old values of components models that determine other components states
 * @param scope
 * @param affectedComponents - list of affected component's models
 */
HeatmapController.prototype.whenFetchHeatmapData =  function (newValues, oldValues, scope, affectedComponents) {

    var preprocessModel = affectedComponents[0],
        runModel = affectedComponents[1],
        fetchBtnFlag = newValues[0],
        fetchSamples = newValues[1],
        fetchSubsets = newValues[2];


    var _toggleRankCriteria = function (runAnalysisModel, noOfSamples, noOfSubsets) {

        runAnalysisModel.subsets = noOfSubsets;

        if (noOfSubsets > 1) {
            runAnalysisModel.params.ranking = 'bval';
        } else {
            runAnalysisModel.params.ranking = 'coef';
        }
        if (noOfSamples === 1) {
            runAnalysisModel.params.ranking = 'mean';
        }
    };


    if (fetchBtnFlag) { // empty preprocess & run result when fetching
        preprocessModel.scriptResults = {};
        runModel.scriptResults = {};
    } else {
        if (fetchSamples > 0 && fetchSubsets > 0) {
            _toggleRankCriteria(runModel, fetchSamples, fetchSubsets);
        }
    }

    preprocessModel.btn.disabled = fetchBtnFlag;
    preprocessModel.btn.disabled = runModel.btn.disabled ? true : fetchSamples <= 1;
    runModel.btn.disabled = fetchBtnFlag;
};

/**
 * Collective behaviour when preprocess data for heatmap
 * @param newValues
 * @param oldValues
 * @param scope
 * @param affectedComponents
 */
HeatmapController.prototype.whenPreprocessHeatmapData = function (newValues, oldValues, scope, affectedComponents) {
    var preprocessBtnFlag = newValues[0],
        fetchModel = affectedComponents[0],
        runModel = affectedComponents[1];

    // toggle fetch & run buttons
    fetchModel.btn.disabled = preprocessBtnFlag;
    runModel.btn.disabled = preprocessBtnFlag;

    if (preprocessBtnFlag) {
        runModel.scriptResults = {};
    }
};

/**
 * Collective behaviour when running heatmap
 * @param newValues
 * @param oldValues
 * @param scope
 * @param affectedComponents
 */
HeatmapController.prototype.whenRunHeatmapAnalysis = function (newValues, oldValues, scope, affectedComponents) {
    var runBtnFlag = newValues[0],
        runResults = newValues[1],
        fetchModel = affectedComponents[0],
        preprocessModel = affectedComponents[1],
        downloadBtns = affectedComponents[2];

    // toggle downloads related buttons
    downloadBtns.disabled = Object.keys(runResults).length === 0;

    // toggle fetch & preprocess buttons
    fetchModel.btn.disabled = runBtnFlag;
    preprocessModel.btn.disabled = runBtnFlag;
};

window.smartRApp.controller('HeatmapController',
    ['$scope', 'smartRUtils', 'commonWorkflowService', HeatmapController]);
