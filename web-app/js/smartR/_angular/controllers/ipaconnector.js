//# sourceURL=ipaconnector.js

'use strict';

window.smartRApp.controller('IpaconnectorController', [
    '$scope',
    'commonWorkflowService',
    'smartRUtils',
    function($scope, commonWorkflowService, smartRUtils) {
        commonWorkflowService.initializeWorkflow('ipaconnector', $scope);

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
            },
            selectedBiomarkers: [],
            scriptResults: {}
        };

        // ------------------------------------------------------------- //
        // IPA                                                           //
        // ------------------------------------------------------------- //
        function clearAll() {
            $scope.differentiallyExpressed = null
        }
        $scope.clearAll = clearAll;

        $scope.$watch('runAnalysis.params', function(newValue) {
            clearAll();
        }, true);

        $scope.runAnalysis = {
            disabled: false,
            running: false,
            params: {
                significanceMeasure: 'pval',
                significanceCutoff: 0.05,
                fcCutoff: 1.5,
            },
            ipaCredentials: {
                username: null,
                password: null,
            },
            ipaConnectionIsSecure: window.location.protocol == 'https:' ? true : false,
            scriptResults: {}
        };

        function listsToSimpleTable(lists) {
            if (lists == null)
                return null;

            var simpleTable;
            if ('header' in lists) {
                simpleTable = {
                    header: lists.header,
                    rows: [],
                };
                if (lists.columnTypes != null) {
                    simpleTable.columnTypes = lists.columnTypes;
                }
                // ensure all columns exist
                for (var c=0; c < simpleTable.header.length; c++) {
                    if (! c in lists)
                        throw "column "+c+" not in lists";
                }
            } else {
                // NOTE this is not ordered
                simpleTable = {
                    header: [],
                    rows: [],
                };
                for (var c in lists) {
                    simpleTable.header.push(c);
                }
            }

            // empty
            if (simpleTable.header.length == 0 || lists[simpleTable.header[0]].length == 0)
                return null;

            for (var i=0; i < lists[simpleTable.header[0]].length; i++) {
                var row = []
                for (var c=0; c < simpleTable.header.length; c++) {
                    row.push(lists[simpleTable.header[c]][i]);
                }
                simpleTable.rows.push(row);
            }

            return simpleTable;
        }

        $scope.$watch('fetch.conceptBoxes.highDimensional', function(newValue) {
            $scope.differentiallyExpressed = null;
        }, true);

        $scope.$watch('runAnalysis.scriptResults', function(newValue) {
            $scope.differentiallyExpressed = listsToSimpleTable($scope.runAnalysis.scriptResults);
        }, true);

        $scope.common = {
            totalSamples: 0,
            numberOfRows: 0,
            subsets: 0
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
