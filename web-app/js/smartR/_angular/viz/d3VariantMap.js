//# sourceURL=d3VariantMap.js

'use strict';

window.smartRApp.directive('variantMap', [
    'smartRUtils',
    'rServeService',
    function(smartRUtils, rServeService) {

        return {
            restrict: 'E',
            scope: {
                data: '=',
                width: '@',
                height: '@'
            },
            link: function (scope, element) {

                scope.$watch('data', function() {
                    $(element[0]).empty();
                    if (! $.isEmptyObject(scope.data)) {
                        smartRUtils.prepareWindowSize(scope.width, scope.height);
                        createVariantMap(scope, element[0]);
                    }
                });
            }
        };

        function createVariantMap(scope, root) {
        }

    }]);
