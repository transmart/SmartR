//# sourceURL=correlation.js

'use strict';

window.smartRApp.controller('CorrelationController',
    ['$scope', 'smartRUtils', 'commonWorkflowService', function($scope, smartRUtils, commonWorkflowService) {

        commonWorkflowService.initializeWorkflow('correlation', $scope);

        $scope.fetch = {
            valid: false,
            disabled: false,
            running: false,
            loaded: false,
            conceptBoxes: {
                datapoints: {concepts: [], valid: false},
                annotations: {concepts: [], valid: true}
            }
        };

        $scope.runAnalysis = {
            disabled: true,
            running: false,
            scriptResults: {},
            params: {
                method: 'pearson'
            }
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

