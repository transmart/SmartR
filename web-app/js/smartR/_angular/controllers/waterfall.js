//# sourceURL=waterfall.js

'use strict';

window.smartRApp.controller('WaterfallController', [
    '$scope',
    'commonWorkflowService',
    'smartRUtils',
    function($scope, commonWorkflowService, smartRUtils) {
        commonWorkflowService.initializeWorkflow('waterfall', $scope);

        $scope.debug = true;

        // ------------------------------------------------------------- //
        // Fetch data                                                    //
        // ------------------------------------------------------------- //
        $scope.fetch = {
            disabled: false,
            running: false,
            loaded: false,
            conceptBoxes: {
                numData: {concepts: [], valid: false},
            }
        };

        $scope.runAnalysis = {
            disabled: true,
            running: false,
            params: {
                selLowRange: '&lt;',
                txtLowRange: 0,
                selHighRange: '&gt;',
                txtHighRange: 0,
            },
            scriptResults: {}
        };


        $scope.$watchGroup(['fetch.running', 'fetch.loaded', 'runAnalysis.running'], function(newValues) {
            var fetchRunning = newValues[0],
                fetchLoaded = newValues[1],
                runAnalysisRunning = newValues[2];

            // clear old results
            if (fetchRunning) {
                $scope.runAnalysis.scriptResults = {};
                $scope.common.subsets = smartRUtils.countCohorts();
            }

            // disable tabs when certain criteria are met
            $scope.fetch.disabled = runAnalysisRunning;
            $scope.runAnalysis.disabled = fetchRunning || !fetchLoaded;
        });
    }]);
