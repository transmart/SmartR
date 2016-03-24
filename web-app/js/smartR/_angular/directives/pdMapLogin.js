//# sourceURL=pdMapLogin.js

'use strict';

window.smartRApp.directive('pdMapLogin', ['$rootScope', function($rootScope) {
    return {
        restrict: 'E',
        scope: {
            criteria: '=',
            subsets : '=',
            login: '=',
            password: '='
        },
        templateUrl: $rootScope.smartRPath +  '/js/smartR/_angular/templates/pdMapLogin.html'
    }
}]);