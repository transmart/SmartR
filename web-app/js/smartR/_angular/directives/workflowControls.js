//# sourceURL=workflowControls.js

'use strict';

window.smartRApp.directive('workflowControls', [
    '$rootScope',
    'smartRUtils',
    function($rootScope, smartRUtils) {
        return {
            restrict: 'E',
            transclude: true,
            templateUrl: $rootScope.smartRPath + '/js/smartR/_angular/templates/workflowControls.html',
            link: function(scope, element) {
                var controls = element.children()[0];
                var scrollbarWidth = smartRUtils.getScrollBarWidth();
                controls.style.bottom = scrollbarWidth + 'px';
                controls.style.right = scrollbarWidth + 105 + 'px';
            }
        };
    }
]);
