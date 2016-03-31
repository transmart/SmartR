//# sourceURL=runButton.js

'use strict';

window.smartRApp.directive('runButton',
    ['$rootScope', 'rServeService', function($rootScope, rServeService) {
        return {
            restrict: 'E',
            scope: {
                disabled: '=',
                storage: '=storeResultsIn',
                script: '@scriptToRun',
                name: '@buttonName',
                serialized: '=',
                arguments: '=argumentsToUse'
            },
            templateUrl: $rootScope.smartRPath +  '/js/smartR/_angular/templates/runButton.html',
            link: function(scope, element) {

                var template_btn = element.children()[0],
                    template_msg = element.children()[1],
                    serialized = scope.serialized;

                template_btn.disabled = scope.disabled;

                scope.$watch('disabled', function (newValue) {
                    template_btn.disabled = newValue;
                }, true);

                var _successCreatePlot = function (response) {
                    template_msg.innerHTML = ''; // empty template
                    if (serialized) { // when results are serialized, we need to deserialized them by
                        // downloading the results files.
                        rServeService.downloadJsonFile(response.executionId, 'heatmap.json').then(
                            function (d) {
                                scope.storage = d.data;
                                scope.disabled = false;
                            }
                        );
                    } else { // results
                        scope.storage = JSON.parse(response.result.artifacts.value);
                        scope.disabled = false;
                    }
                };

                var _failCreatePlot = function (response) {
                    template_msg.style.color = 'red';
                    template_msg.innerHTML = 'Failure: ' + response.statusText;
                    scope.disabled = false;
                };

                template_btn.onclick = function() {

                    scope.storage = {};
                    scope.disabled = true;
                    template_msg.innerHTML = 'Creating plot, please wait <span class="blink_me">_</span>';

                    rServeService.startScriptExecution({
                        taskType: scope.script,
                        arguments: scope.arguments
                    }).then(
                        _successCreatePlot,
                        _failCreatePlot
                    );
                };
            }
        };
    }]);
