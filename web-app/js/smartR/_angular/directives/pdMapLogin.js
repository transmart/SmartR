//# sourceURL=pdMapLogin.js

window.smartRApp.directive('pdMapLogin', ['$rootScope', function($rootScope) {
    return {
        restrict: 'E',
        scope: {
            criteria: '=',
            subsets : '=',
            login: '=',
            password: '=',
            pdMapLinkEnabled: '='
        },
        templateUrl: $rootScope.smartRPath +  '/js/smartR/_angular/templates/pdMapLogin.html'
    }
}]);