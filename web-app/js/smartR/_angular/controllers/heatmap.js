//# sourceURL=heatmap.js

window.smartRApp.controller('HeatmapController',
    ['$scope', 'smartRUtils', 'commonWorkflowService', function($scope, smartRUtils, commonWorkflowService) {

        commonWorkflowService.initializeWorkflow('heatmap', $scope);

        // ------------------------------------------------------------- //
        // Fetch data                                                    //
        // ------------------------------------------------------------- //
        $scope.fetch = {
            conceptBoxes : {
                highDimensional : [],
                numerical : [],
                categorical : []
            },
            selectedBiomarkers : [],
            summaryData : {}
        };

        // ------------------------------------------------------------- //
        // Preprocess                                                    //
        // ------------------------------------------------------------- //
        $scope.preprocess = {
            params :  {
                aggregateProbes : false
            },
            summaryData : {}
        };

        // ------------------------------------------------------------- //
        // Run Heatmap                                                   //
        // ------------------------------------------------------------- //
        $scope.runAnalysis = {
            params: {
                max_row : 100,
                sorting : 'nodes',
                ranking : 'coef',
                pdMapLinkEnabled: false,
                pdMapLogin: 'Enter Login',
                pdMapPassword: 'Enter Password'
            },
            subsets : 0,
            scriptResults : {}
        };
    }]);
