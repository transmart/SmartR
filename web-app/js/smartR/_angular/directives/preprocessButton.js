//# sourceURL=preprocessButton.js

window.smartRApp.directive('preprocessButton', ['rServeService',
    function(rServeService) {
    return {
        restrict: 'E',
        scope: {
            disabled: '=',
            params: '=',
            showSummaryStats: '=',
            summaryData: '='
        },
        template:
            '<input type="button" value="Preprocess" class="heim-action-button">' +
            '<span style="padding-left: 10px;"></span>',
        link: function(scope, element) {

            var template_btn = element.children()[0];
            var template_msg = element.children()[1];

            template_btn.disabled = scope.disabled;

            scope.$watch('disabled', function (newValue) {
                template_btn.disabled = newValue;
            }, true);

            scope.$watch('summaryData', function (newSummaryData) {
                if (newSummaryData.hasOwnProperty('allSamples')) {
                    // when everything is retrieved
                    scope.disabled = false;
                }
            }, true);

            template_btn.onclick = function() {

                var _init = function () {
                        scope.summaryData = {}; // reset
                        scope.disabled = true;
                        template_msg.innerHTML = 'Preprocessing, please wait <span class="blink_me">_</span>';
                    },

                    _args = {aggregate:scope.params.aggregateProbes},

                    _preprocessData = function (_args) {
                        return rServeService.preprocess(_args).then(function (msg){
                            return msg;
                        });
                    },

                    _finishedPreprocessed = function (msg) {
                        template_msg.innerHTML = 'Success: ' + msg;
                        scope.disabled = false;
                    },

                    _afterDataPreprocessed = function (msg) {
                        if (!scope.showSummaryStats) {
                            return _finishedPreprocessed(msg);
                        }
                        template_msg.innerHTML = 'Execute summary statistics, please wait <span class="blink_me">_</span>';

                        return  rServeService.executeSummaryStats('preprocess')
                            .then (function (data) {
                                scope.summaryData = data.result;
                                template_msg.innerHTML = 'Success: ' + data.msg;
                            }, function(msg) {
                                template_msg.innerHTML = 'Failure: ' + msg;
                            })
                    };

                _init();

                _preprocessData(_args)
                    .then(_afterDataPreprocessed);

            };
        }
    };
}]);
