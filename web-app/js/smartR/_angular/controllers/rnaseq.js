//# sourceURL=rnaseq.js

'use strict';

window.smartRApp.controller('RNASeqController', [
    '$scope',
    'commonWorkflowService',
    'smartRUtils',
    function($scope, commonWorkflowService, smartRUtils) {
        commonWorkflowService.initializeWorkflow('rnaseq', $scope);

        $scope.debug = true;

        // ------------------------------------------------------------- //
        // Fetch data                                                    //
        // ------------------------------------------------------------- //
        $scope.fetch = {
            disabled: false,
            running: false,
            loaded: false,
            conceptBoxes: {
                highDimensional: {concepts: [], valid: false},
                cat_grouping: {concepts: [], valid: true},
                num_grouping: {concepts: [], valid: true}
            }
        };

        $scope.runAnalysis = {
            disabled: true,
            running: false,
            params: {
                analysis_type: 'unpaired_group'
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
            }

            // disable tabs when certain criteria are met
            $scope.fetch.disabled = runAnalysisRunning;
            $scope.runAnalysis.disabled = fetchRunning || !fetchLoaded;
        });
    }]);
