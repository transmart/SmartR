//# sourceURL=sortingCriteria.js

window.smartRApp.directive('sortingCriteria', ['$rootScope', function($rootScope) {
    return {
        restrict: 'E',
        scope: {
            criteria : '=',
            subsets : '='
        },
        templateUrl: $rootScope.smartRPath +  '/js/smartR/_angular/templates/sortingCriteria.html'
    }
}]);
