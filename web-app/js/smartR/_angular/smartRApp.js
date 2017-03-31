//# sourceURL=smartRApp.js

'use strict';

window.smartRApp = angular.module('smartRApp', ['ngRoute', 'door3.css', 'ipaApi'])
    .config(['$httpProvider', function($httpProvider) {
        //initialize get if not there
        if (!$httpProvider.defaults.headers.get) {
            $httpProvider.defaults.headers.get = {};
        }
        //disable IE ajax request caching
        $httpProvider.defaults.headers.get['If-Modified-Since'] = 'Mon, 26 Jul 1997 05:00:00 GMT';
    }])
    .run(function($rootScope, $http) {
        // get plugin context path and put it in root scope
        $http.get(pageInfo.basePath + '/SmartR/smartRContextPath').then(
            function(d) { $rootScope.smartRPath = d.data; },
            function(msg) { throw 'Error: ' + msg; }
        );
    });
