//# sourceURL=logitregression.js

'use strict';

window.smartRApp.controller('LogitRegressionController', [
    '$scope',
    'smartRUtils',
    'commonWorkflowService',
    function($scope, smartRUtils, commonWorkflowService) {

        commonWorkflowService.initializeWorkflow('logitregression', $scope);

        $scope.fetch = {
            disabled: true,
            running: false,
            loaded: false,
            conceptBoxes: {
                datapoints: {concepts: [], valid: false},
                annotations: {concepts: [], valid: true}
            }
        };

        $scope.runAnalysis = {
            params: {
                transformationx: 'raw',
                transformationy: 'raw'
            },
            disabled: true,
            running: false,
            scriptResults: {}
        };


        // $scope.$watch('runAnalysis.params.transformation', function(newValue, oldValue) {
        //     ""
        // });


        $scope.$watchGroup(['fetch.running', 'runAnalysis.running'], function(newValues) {
            var fetchRunning = newValues[0],
                runAnalysisRunning = newValues[1];

            // clear old results
            if (fetchRunning) {
                $scope.runAnalysis.scriptResults = {};
            }

            // disable tabs when certain criteria are not met
            $scope.fetch.disabled = runAnalysisRunning;
            $scope.runAnalysis.disabled = fetchRunning || !$scope.fetch.loaded;
        });

    }]);

