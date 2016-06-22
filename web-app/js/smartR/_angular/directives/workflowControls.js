//# sourceURL=workflowControls.js

'use strict';

window.smartRApp.directive('workflowControls', [
    '$rootScope',
    function($rootScope) {
        return {
            restrict: 'E',
            transclude: true,
            templateUrl: $rootScope.smartRPath + '/js/smartR/_angular/templates/workflowControls.html'
        };
    }
]);
