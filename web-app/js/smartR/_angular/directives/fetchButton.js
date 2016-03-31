//# sourceURL=fetchButton.js

'use strict';

window.smartRApp.directive('fetchButton',
    ['$rootScope', 'rServeService', 'smartRUtils',
        function($rootScope, rServeService, smartRUtils) {
        return {
            restrict: 'E',
            scope: {
                disabled: '=',
                conceptMap: '=',
                biomarkers: '=?',
                showSummaryStats: '=?',
                summaryData: '=?',
                allSamples: '=?',
                subsets: '=?'
            },
            templateUrl: $rootScope.smartRPath +  '/js/smartR/_angular/templates/fetchButton.html',
            link: function(scope, element) {

                var template_btn = element.children()[0],
                    template_msg = element.children()[1];

                template_btn.disabled = scope.disabled;

                scope.$watch('disabled', function (newValue) {
                    template_btn.disabled = newValue;
                }, true);

                if (angular.isDefined(scope.summaryData)) {
                    scope.$watch('summaryData', function (newSummaryData) {
                        if (newSummaryData.hasOwnProperty('allSamples')) {
                            // when everything is retrieved
                            scope.allSamples = newSummaryData.allSamples;
                            scope.subsets = newSummaryData.subsets;
                            scope.disabled = false;
                        }
                    }, true);
                }

                template_btn.onclick = function() {

                    var _init = function () {
                            scope.summaryData = {}; // reset
                            scope.allSamples = 0;
                            scope.subsets = 0;
                            scope.disabled = true;
                            template_msg.innerHTML = 'Fetching data, please wait <span class="blink_me">_</span>';
                        },

                        _getDataConstraints = function (biomarkers) {
                            if (typeof biomarkers !== 'undefined' && biomarkers.length > 0) {
                                var searchKeywordIds = biomarkers.map(function(biomarker) {
                                    return String(biomarker.id);
                                });
                                return dataConstraints = {
                                    search_keyword_ids: {
                                        keyword_ids: searchKeywordIds
                                    }
                                };
                            }
                        },

                        _fetchDataIntoRSession = function (conceptKeys, dataConstraints) {
                            return rServeService.loadDataIntoSession(conceptKeys, dataConstraints)
                                .then(function(msg) {
                                    return msg;
                                });
                        },

                        _finishedFetching = function (msg) {
                            template_msg.innerHTML = 'Success: ' + msg;
                            scope.disabled = false;
                        },

                        _afterDataFetched = function (msg) {

                            if (!scope.showSummaryStats) {
                                return _finishedFetching(msg);
                            }

                            template_msg.innerHTML =
                                'Execute summary statistics, please wait <span class="blink_me">_</span>';

                            return rServeService.executeSummaryStats('fetch')
                                .then (function (data) {
                                    scope.summaryData = data.result;
                                    template_msg.innerHTML = 'Success: ' + data.msg;
                                }, function(msg) {
                                    template_msg.innerHTML = 'Failure: ' + msg;
                                    scope.disabled = false;
                                })
                        },
                        conceptKeys = smartRUtils.conceptBoxMapToConceptKeys(scope.conceptMap),
                        dataConstraints = _getDataConstraints(scope.biomarkers);

                    _init();

                    _fetchDataIntoRSession(conceptKeys, dataConstraints)
                        .then(_afterDataFetched);

                }; // end onclick
            }
        };
    }]);
