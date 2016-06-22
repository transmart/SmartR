//# sourceURL=heatmap.js

'use strict';

window.smartRApp.controller('HeatmapController', [
    '$scope',
    'commonWorkflowService',
    'smartRUtils',
    function($scope, commonWorkflowService, smartRUtils) {

        commonWorkflowService.initializeWorkflow('heatmap', $scope);

        // ------------------------------------------------------------- //
        // Fetch data                                                    //
        // ------------------------------------------------------------- //
        $scope.fetch = {
            disabled: false,
            running: false,
            loaded: false,
            conceptBoxes: {
                highDimensional: {concepts: [], valid: false},
                numeric: {concepts: [], valid: true},
                categoric: {concepts: [], valid: true}
            },
            selectedBiomarkers: [],
            scriptResults: {}
        };

        // ------------------------------------------------------------- //
        // Preprocess                                                    //
        // ------------------------------------------------------------- //
        $scope.preprocess = {
            disabled: true,
            running: false,
            params:  {
                aggregate: false
            },
            scriptResults: {}
        };

        // ------------------------------------------------------------- //
        // Run Heatmap                                                   //
        // ------------------------------------------------------------- //
        $scope.runAnalysis = {
            disabled: true,
            running: false,
            params: {
                max_row: 100,
                sorting: 'nodes',
                ranking: 'coef',
                geneCardsAllowed: false,
            },
            download: {
                disabled: true
            },
            scriptResults: {}
        };

        $scope.common = {
            totalSamples: 0,
            numberOfRows: 0,
            subsets: 0
        };

        $scope.$watchGroup(['fetch.running', 'preprocess.running', 'runAnalysis.running'], function(newValues) {
            var fetchRunning = newValues[0],
                preprocessRunning = newValues[1],
                runAnalysisRunning = newValues[2];

            // clear old results
            if (fetchRunning) {
                $scope.preprocess.scriptResults = {};
                $scope.runAnalysis.scriptResults = {};
                $scope.runAnalysis.params.ranking = '';
                $scope.common.subsets = smartRUtils.countCohorts();
            }

            // clear old results
            if (preprocessRunning) {
                $scope.runAnalysis.scriptResults = {};
                $scope.runAnalysis.params.ranking = '';
                $scope.common.subsets = smartRUtils.countCohorts();
            }

            // disable tabs when certain criteria are not met
            $scope.fetch.disabled = preprocessRunning || runAnalysisRunning;
            $scope.preprocess.disabled = fetchRunning || runAnalysisRunning || !$scope.fetch.loaded;
            $scope.runAnalysis.disabled = fetchRunning || preprocessRunning || !$scope.fetch.loaded;

            // disable buttons when certain criteria are not met
            $scope.runAnalysis.download.disabled = runAnalysisRunning ||
                $.isEmptyObject($scope.runAnalysis.scriptResults);

            // set ranking criteria
            if (!fetchRunning &&
                !preprocessRunning &&
                $scope.common.totalSamples < 2 &&
                $scope.runAnalysis.params.ranking === '') {
                $scope.runAnalysis.params.ranking = 'mean';
            } else if (!fetchRunning &&
                       !preprocessRunning &&
                       $scope.common.subsets < 2 &&
                       $scope.runAnalysis.params.ranking === '') {
                $scope.runAnalysis.params.ranking = 'coef';
            } else if (!fetchRunning &&
                       !preprocessRunning &&
                       $scope.common.subsets > 1 &&
                       $scope.runAnalysis.params.ranking === '') {
                $scope.runAnalysis.params.ranking = 'adjpval';
            }
        });
    }]);
