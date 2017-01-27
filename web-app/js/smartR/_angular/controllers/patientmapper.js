//# sourceURL=patientmapper.js

'use strict';

window.smartRApp.controller('PatientmapperController',
    ['$scope', 'smartRUtils', 'commonWorkflowService', function($scope, smartRUtils, commonWorkflowService) {

        commonWorkflowService.initializeWorkflow('patientmapper', $scope);

        $scope.fetch = {
            disabled: false,
            running: false,
            loaded: false,
            conceptBoxes: {
                source: {concepts: [], valid: false},
                target: {concepts: [], valid: false}
            }
        };

        $scope.runAnalysis = {
            disabled: true,
            running: false,
            scriptResults: {},
            params: {}
        };


        $scope.$watch('runAnalysis.scriptResults', function(results) {
            if (!$.isEmptyObject(results)) {
                clearQuery();
                results.cohortNodes.forEach(function(d) {
                    var panel = document.querySelector('.panelModel[subset="' + d.subset + '"] .panelBoxList');
                    var concept = {
                        name: d.name,
                        key: d.key,
                        level: d.level,
                        tooltip: d.tooltip,
                        dimcode: d.fullname,
                        value: {
                            mode: 'novalue',
                            operator: 'LT',
                            highlowselect: 'N',
                        },
                        oktousevalues: 'N',
                        nodeType: 'alphaicon',
                        visualattributes: d.visualAttributes,
                        applied_path: '@',
                        modifiedNode: {}
                    };
                    selectConcept(createPanelItemNew(panel, concept));
                });
                appendQueryPanelInto(1);
                appendQueryPanelInto(2);
                runAllQueries(function() {
                    $.ajax({
                        url: pageInfo.basePath + '/chart/clearGrid',
                        method: 'POST',
                        data: {
                            charttype: 'cleargrid',
                        }
                    });
                });
            }
        });
    }]);

