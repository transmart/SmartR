
'use strict';

window.smartRApp.controller('BoxplotController', [
    '$scope',
    'smartRUtils',
    'commonWorkflowService',
    function($scope, smartRUtils, commonWorkflowService) {

        commonWorkflowService.initializeWorkflow('boxplot', $scope);

        $scope.fetch = {
            running: false,
            disabled: false,
            loaded: false,
            conceptBoxes: {
                datapoints: {concepts: [], valid: false}
            }
        };

        $scope.runAnalysis = {
            running: false,
            disabled: true,
            scriptResults: {},
            params: {}
        };

        $scope.$watchGroup(['fetch.running', 'runAnalysis.running'],
            function(newValues) {
                var fetchRunning = newValues[0],
                    runAnalysisRunning = newValues[1];

                // clear old results
                if (fetchRunning) {
                    $scope.runAnalysis.scriptResults = {};
                }

                // disable tabs when certain criteria are not met
                $scope.fetch.disabled = runAnalysisRunning;
                $scope.runAnalysis.disabled = fetchRunning || !$scope.fetch.loaded;
            }
        );

    }]);

