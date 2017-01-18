//# sourceURL=fetchButton.js

'use strict';

window.smartRApp.directive('fetchButton', [
    '$rootScope',
    'rServeService',
    'smartRUtils',
    function($rootScope, rServeService, smartRUtils) {
        return {
            restrict: 'E',
            scope: {
                conceptMap: '=',
                loaded: '=?',
                running: '=?',
                biomarkers: '=?',
                showSummaryStats: '=?',
                summaryData: '=?',
                allSamples: '=?',
                numberOfRows: '=?',
                allowedCohorts: '=',
                projection: '@?',
                disabled: '=?',
                message: '=?',
                hasPreprocessTab: '=?',
                callback: '&?',
            },
            templateUrl: $rootScope.smartRPath +  '/js/smartR/_angular/templates/fetchButton.html',
            link: function(scope, element) {
                var template_btn = element.children()[0],
                    template_msg = element.children()[1];

                var _onSuccess = function() {
                    if (scope.hasPreprocessTab) {
                        scope.message = 'Task complete! Go to the "Preprocess" or "Run Analysis" tab to continue.';
                    } else {
                        scope.message = 'Task complete! Go to the "Run Analysis" tab to continue.';
                    }
                    scope.loaded = true;
                    scope.disabled = false;
                    scope.running = false;
                    scope.$apply();
                };

                var _onFailure = function(msg) {
                    scope.message = 'Error: ' + msg;
                    scope.loaded = false;
                    scope.disabled = false;
                    scope.running = false;
                    scope.$apply();
                };

                // we add this conditional $watch because there is some crazy promise resolving for allSamples
                // going on. This is a workaround which observes allSamples and uses it as criteria for successful
                // completion. FIXME
                scope.$watch('summaryData', function(newValue) {
                    if (scope.summaryData &&
                        scope.showSummaryStats &&
                        scope.running &&
                        Object.keys(newValue).indexOf('subsets') !== -1) {
                        scope.allSamples = newValue.allSamples;
                        scope.numberOfRows = newValue.numberOfRows;
                        _onSuccess();
                    }
                }, true);

                scope.$watch('message', function(newValue) {
                    if (newValue) {
                        template_msg.innerHTML = newValue;
                    } else {
                        template_msg.innerHTML = '';
                    }
                });

                var _getDataConstraints = function (biomarkers) {
                    if (typeof biomarkers !== 'undefined' && biomarkers.length > 0) {
                        var searchKeywordIds = biomarkers.map(function(biomarker) {
                            return String(biomarker.id);
                        });
                        return {
                            search_keyword_ids: {
                                keyword_ids: searchKeywordIds
                            }
                        };
                    }
                };

                var _showSummaryStats = function() {
                    scope.message = 'Executing summary statistics, please wait <span class="blink_me">_</span>';
                    rServeService.executeSummaryStats('fetch')
                        .then(
                            function(data) { scope.summaryData = data.result; }, // this will trigger $watch
                            _onFailure
                        );
                };

                template_btn.onclick = function() {
                    scope.disabled = true;
                    scope.message = 'Fetching data, please wait <span class="blink_me">_</span>';

                    scope.summaryData = {};
                    scope.allSamples = 0;
                    scope.loaded = false;
                    scope.running = true;
                    var deleteReq = rServeService.deleteSessionFiles(); // cleanup our working directory
                    var cohorts = smartRUtils.countCohorts();

                    if (cohorts === 0) {
                        _onFailure('No cohorts selected!');
                        return;
                    }

                    if (scope.allowedCohorts.indexOf(cohorts) === -1) {
                        _onFailure('This workflow requires ' + scope.allowedCohorts +
                            ' cohort(s), but you selected ' + cohorts);
                        return;
                    }

                    for (var conceptGroup in scope.conceptMap) {
                        if (scope.conceptMap.hasOwnProperty(conceptGroup) && !scope.conceptMap[conceptGroup].valid) {
                            _onFailure('Your data do not match the requirements! All fields must be green.');
                            return;
                        }
                    }

                    var conceptKeys = smartRUtils.conceptBoxMapToConceptKeys(scope.conceptMap);
                    if ($.isEmptyObject(conceptKeys) && !scope.biomarkers.length) {
                        _onFailure('No concepts or genes selected!');
                        return;
                    }

                    var dataConstraints = _getDataConstraints(scope.biomarkers);
                    deleteReq.then(
                        function() {
                            if ($.isEmptyObject(conceptKeys)) {
                                if (scope.callback) {
                                    scope.callback()().then(
                                        scope.showSummaryStats ? _showSummaryStats : _onSuccess,
                                        _onFailure
                                    );
                                } else {
                                    _onSuccess();
                                }
                            } else {
                                rServeService.loadDataIntoSession(conceptKeys, dataConstraints, scope.projection).then(
                                    function() {
                                        if (scope.callback) {
                                            scope.callback()().then(
                                                scope.showSummaryStats ? _showSummaryStats : _onSuccess,
                                                _onFailure
                                            );
                                        } else {
                                            _onSuccess();
                                        }
                                    },
                                    _onFailure
                                );
                            }
                        }, _onFailure
                    );
                };
            }
        };
    }]);
