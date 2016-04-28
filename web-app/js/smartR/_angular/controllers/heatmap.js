//# sourceURL=heatmap.js

'use strict';

window.smartRApp.controller('HeatmapController', [
    '$scope',
    'commonWorkflowService',
    function($scope, commonWorkflowService) {

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
                numerical: {concepts: [], valid: true},
                categorical: {concepts: [], valid: true}
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
            subsets: 0,
            totalSamples: 0,
            numberOfRows: 0
        };

        $scope.$watchGroup(['fetch.running', 'preprocess.running', 'runAnalysis.running'],
            function(newValues) {
                var fetchRunning = newValues[0],
                    preprocessRunning = newValues[1],
                    runAnalysisRunning = newValues[2];

                // clear old results
                if (fetchRunning) {
                    $scope.preprocess.scriptResults = {};
                    $scope.runAnalysis.scriptResults = {};
                }

                // clear old results
                if (preprocessRunning) {
                    $scope.runAnalysis.scriptResults = {};
                }

                // disable tabs when certain criteria are not met
                $scope.fetch.disabled = preprocessRunning || runAnalysisRunning;
                $scope.preprocess.disabled = fetchRunning || runAnalysisRunning || !$scope.fetch.loaded ||
                    $scope.common.totalSamples <= 1;
                $scope.runAnalysis.disabled = fetchRunning || preprocessRunning || !$scope.fetch.loaded;

                // disable buttons when certain criteria are not met
                $scope.runAnalysis.download.disabled = runAnalysisRunning ||
                    $.isEmptyObject($scope.runAnalysis.scriptResults);

                // set ranking criteria
                if ($scope.common.totalSamples < 2 &&
                        $scope.runAnalysis.params.ranking !== 'median') {
                    $scope.runAnalysis.params.ranking = 'mean';
                } else if ($scope.common.subsets < 2 &&
                        $scope.runAnalysis.params.ranking !== 'variance' &&
                        $scope.runAnalysis.params.ranking !== 'range') {
                    $scope.runAnalysis.params.ranking = 'coef';
                }
            }
        );
    }]);
