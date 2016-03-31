//# sourceURL=downloadResultsButton.js

window.smartRApp.directive('downloadResultsButton', ['rServeService', function(rServeService)
{
    function downloadFile(data) {
        var link = jQuery('<a/>')
            .attr('href', rServeService.urlForFile(data.executionId, 'analysis_data.zip'))
            .attr('download', 'analysis_data.zip')
            .css('display', 'none');
        jQuery('body').append(link);
        link[0].click();
        link.remove();
    }

    return {
        restrict: 'E',
        scope: {
            disabled: '='
        },
        template:
            '<input type="button" value="Download" class="heim-action-button">' +
            '<span style="padding-left: 10px;"></span>',
        link: function(scope, element) {

            var template_btn = element.children()[0];
            var template_msg = element.children()[1];

            template_btn.disabled = scope.disabled;

            scope.$watch('disabled', function (newValue) {
                template_btn.disabled = newValue;
            }, true);

            template_btn.onclick = function() {

                template_msg.innerHTML = 'Download data, please wait <span class="blink_me">_</span>';

                rServeService.startScriptExecution({
                    taskType: 'downloadData',
                    arguments: {}
                }).then(
                    function (data){
                        // download file
                        template_msg.innerHTML = '';
                        downloadFile(data);
                    },
                    function (msg){
                        template_msg.innerHTML = 'Failure: ' + msg;
                    }
                )
            };
        }
    };
}]);
