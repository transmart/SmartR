//# sourceURL=pdMapLogin.js

window.smartRApp.directive('pdMapLogin', ['$rootScope', function($rootScope) {
    return {
        restrict: 'E',
        scope: {
            criteria: '=',
            subsets : '=',
            login: '=',
            password: '=',
            enabled: '='
        },
        templateUrl: $rootScope.smartRPath +  '/js/smartR/_angular/templates/pdMapLogin.html',
        link: function(scope, element) {
            scope.$watchGroup(
                ['enabled', 'login', 'password', 'criteria', 'subsets'], function(newValues) {
                    var enabled = JSON.parse(newValues[0]);
                    var login = newValues[1];
                    var password = newValues[2];
                    var criteria = newValues[3];
                    var subsets = parseInt(newValues[4]);

                    if (subsets > 1 && ['ttest', 'pval', 'adjpval', 'bval', 'logfold'].indexOf(criteria) !== -1) {
                        element[0].querySelector('#sr-pdMapSelect').disabled = false;
                    } else {
                        element[0].querySelector('#sr-pdMapSelect').disabled = true;
                        scope.enabled = false;
                    }

                    if (enabled) {
                        element[0].querySelector('#sr-pdmap-login').disabled = false;
                        element[0].querySelector('#sr-pdmap-password').disabled = false;
                    } else {
                        element[0].querySelector('#sr-pdmap-login').disabled = true;
                        element[0].querySelector('#sr-pdmap-password').disabled = true;
                    }
            });
        }
    }
}]);