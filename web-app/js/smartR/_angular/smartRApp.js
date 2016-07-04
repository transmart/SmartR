//# sourceURL=smartRApp.js

'use strict';

/**
 * Extends Array to be able to get min value
 * @returns {*}
 */
Array.prototype.min = function() {
    var min = this[0];

    for ( var i = 0; i < this.length; i++ )
        if ( this[i] < min ) min = this[i];

    return min;
};

/**
 * Extends Array to be able to get max value
 * @returns {*}
 */
Array.prototype.max = function() {
    var max = this[0];

    for ( var i = 0; i < this.length; i++ )
        if ( this[i] > max ) max = this[i];

    return max;
};


window.smartRApp = angular.module('smartRApp', ['ngRoute', 'door3.css'])
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
