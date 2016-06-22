//# sourceURL=runButton.js

'use strict';

window.smartRApp.directive('runButton', [
    '$rootScope',
    'rServeService',
    'smartRUtils',
    function($rootScope, rServeService, smartRUtils) {
        return {
            restrict: 'E',
            scope: {
                running: '=?',
                storage: '=storeResultsIn',
                script: '@scriptToRun',
                name: '@buttonName',
                filename: '@?',
                params: '=?argumentsToUse'
            },
            templateUrl: $rootScope.smartRPath + '/js/smartR/_angular/templates/runButton.html',
            link: function(scope, element) {
                var params = scope.params ? scope.params : {};

                var template_btn = element.children()[0],
                    template_msg = element.children()[1];

                var _onSuccess = function(data) {
                    scope.storage = data;
                    template_msg.innerHTML = '';
                    template_btn.disabled = false;
                    scope.disabled = false;
                    scope.running = false;
                };

                var _onFail = function(msg) {
                    template_msg.innerHTML = 'Error: ' + msg;
                    template_btn.disabled = false;
                    scope.disabled = false;
                    scope.running = false;
                };

                var _prepareResults = function(response) {
                    if (scope.filename) {
                        // when filename is specified it is assumed that results are serialized
                        // if results are serialized, we need to deserialized them by
                        // downloading the results files.
                        rServeService.downloadJsonFile(response.executionId, scope.filename).then(
                            function(d) { _onSuccess(d.data); },
                            _onFail
                        );
                    } else {
                        _onSuccess(JSON.parse(response.result.artifacts.value));
                    }
                };

                template_btn.onclick = function() {
                    smartRUtils.cleanUp();
                    template_btn.disabled = true;
                    scope.storage = {};
                    scope.disabled = true;
                    scope.running = true;
                    template_msg.innerHTML = 'Creating plot, please wait <span class="blink_me">_</span>';

                    rServeService.startScriptExecution({
                        taskType: scope.script,
                        arguments: params
                    }).then(
                        _prepareResults,
                        _onFail
                    );
                };
            }
        };
    }]);
