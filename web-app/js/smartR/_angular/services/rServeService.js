//# sourceURL=rServeService.js

'use strict';

window.smartRApp.factory('rServeService', [
    'smartRUtils',
    '$q',
    '$http',
    function(smartRUtils, $q, $http) {

        var service = {};

        var NOOP_ABORT = function() {};
        var TIMEOUT = 10000; // 10 s
        var CHECK_DELAY = 500; // 0.5 s
        var SESSION_TOUCH_DELAY = 60000; // 1 min

        /* we only support one session at a time */

        var state = {
            currentRequestAbort: NOOP_ABORT,
            sessionId: null,
            touchTimeout: null // for current session id
        };

        var workflow = '';
        /* returns a promise with the session id and
         * saves the session id for future calls */
        service.startSession = function(name) {
            workflow = name;
            var request = $http({
                url: pageInfo.basePath + '/RSession/create',
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json'
                },
                config: {
                    timeout: TIMEOUT
                },
                data: {
                    workflow: workflow
                }
            });

            return $q(function(resolve, reject) {
                request.then(
                    function(response) {
                        state.sessionId = response.data.sessionId;
                        rServeService_scheduleTouch();
                        resolve();
                    },
                    function(response) { reject(response.statusText);}
                );
            });
        };

        service.touch = function(sessionId) {
            if (sessionId !== state.sessionId) {
                return;
            }

            var touchRequest = $http({
                url: pageInfo.basePath + '/RSession/touch',
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json'
                },
                config: {
                    timeout: TIMEOUT
                },
                data: {
                    sessionId: sessionId
                }
            });

            touchRequest.finally(function() {
                rServeService_scheduleTouch(); // schedule another
            });
        };

        function rServeService_scheduleTouch() {
            window.clearTimeout(state.touchTimeout);
            state.touchTimeout = window.setTimeout(function() {
                service.touch(state.sessionId);
            }, SESSION_TOUCH_DELAY);
        }

        service.deleteSessionFiles = function(sessionId) {
            sessionId = sessionId || state.sessionId;

            return $http({
                url: pageInfo.basePath + '/RSession/deleteFiles',
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json'
                },
                config: {
                    timeout: TIMEOUT
                },
                data: {
                    sessionId: sessionId
                }
            });
        };

        service.destroySession = function(sessionId) {
            sessionId = sessionId || state.sessionId;

            if (!sessionId) {
                return;
            }

            var request = $http({
                url: pageInfo.basePath + '/RSession/delete',
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json'
                },
                config: {
                    timeout: TIMEOUT
                },
                data: {
                    sessionId: sessionId
                }
            });

            return request.finally(function() {
                if (state.sessionId === sessionId) {
                    service.abandonCurrentSession();
                }
            });
        };

        service.abandonCurrentSession = function() {
            window.clearTimeout(state.touchTimeout);
            state.sessionId = null;
        };

        service.destroyAndStartSession = function(workflowName) {
            $q.when(service.destroySession()).then(function() {
                service.startSession(workflowName);
            });
        };

        /*
         * taskData = {
         *     arguments: { ... },
         *     taskType: 'fetchData' or name of R script minus .R,
         *     phase: 'fetch' | 'preprocess' | 'run',
         * }
         */
        service.startScriptExecution = function(taskDataOrig) {

            var taskData = $.extend({}, taskDataOrig); // clone the thing
            state.currentRequestAbort();

            var canceler = $q.defer();
            var runRequest = $http({
                url: pageInfo.basePath + '/ScriptExecution/run',
                method: 'POST',
                timeout: canceler.promise,
                responseType: 'json',
                data: {
                    sessionId: state.sessionId,
                    arguments: taskData.arguments,
                    taskType: taskData.taskType,
                    workflow: workflow
                }
            });

            runRequest.finally(function() {
                state.currentRequestAbort = NOOP_ABORT;
            });

            state.currentRequestAbort = function() { canceler.resolve(); };

            /* schedule checks */
            var promise = $q(function(resolve, reject) {
                runRequest.then(
                    function(response) {
                        if (!response.data) {
                            console.error('Unexpected response:', response);
                        }
                        taskData.executionId = response.data.executionId;
                        _checkStatus(taskData.executionId, resolve, reject);
                    },
                    function(response) {
                        reject(response.statusText);
                    }
                );
            });

            promise.cancel = function() {
                // calling this method should by itself resolve the promise
                state.currentRequestAbort();
            };

            // no touching necessary when a task is running
            window.clearTimeout(state.touchTimeout);
            promise.finally(rServeService_scheduleTouch.bind(this));

            return promise;
        };

        /* aux function of _startScriptExecution. Needs to follow its contract
         * with respect to the fail and success result of the promise */
        function _checkStatus(executionId, resolve, reject) {
            var canceler = $q.defer();
            var statusRequest = $http({
                method: 'GET',
                timeout: canceler.promise,
                url: pageInfo.basePath + '/ScriptExecution/status' +
                    '?sessionId=' + state.sessionId +
                    '&executionId=' + executionId
            });

            statusRequest.finally(function() { state.currentRequestAbort = NOOP_ABORT; });
            state.currentRequestAbort = function() { canceler.resolve(); };

            statusRequest.then(
                function (d) {
                    if (d.data.state === 'FINISHED') {
                        d.data.executionId = executionId;
                        resolve(d.data);
                    } else if (d.data.state === 'FAILED') {
                        reject(d.data.result.exception);
                    } else {
                        // else still pending
                        window.setTimeout(function() {
                            _checkStatus(executionId, resolve, reject);
                        }, CHECK_DELAY);
                    }
                },
                function(response) { reject(response.statusText); }
            );
        }

        service.downloadJsonFile = function(executionId, filename) {
            return $http({
                method: 'GET',
                url: this.urlForFile(executionId, filename)
            });
        };


        service.urlForFile = function(executionId, filename) {
            return pageInfo.basePath +
                '/ScriptExecution/downloadFile?sessionId=' + state.sessionId +
                '&executionId=' + executionId + '&filename=' + filename;
        };

        service.loadDataIntoSession = function(conceptKeys, dataConstraints, projection) {
            projection = typeof projection === 'undefined' ? 'log_intensity' : projection; // default to log_intensity
            return $q(function(resolve, reject) {
                smartRUtils.getSubsetIds().then(
                    function(subsets) {
                        var _arg = {
                            conceptKeys: conceptKeys,
                            resultInstanceIds: subsets,
                            projection: projection
                        };

                        if (typeof dataConstraints !== 'undefined') {
                            _arg.dataConstraints = dataConstraints;
                        }

                        service.startScriptExecution({
                            taskType: 'fetchData',
                            arguments: _arg
                        }).then(
                            resolve,
                            function(response) { reject(response); }
                        );
                    },
                    function() {
                        reject('Could not create subsets!');
                    }
                );
            });
        };

        service.executeSummaryStats = function(phase, projection) {
            projection = typeof projection === 'undefined' ? 'log_intensity' : projection; // default to log_intensity
            return $q(function(resolve, reject) {
                service.startScriptExecution({
                    taskType: 'summary',
                    arguments: {
                        phase: phase,
                        projection: projection // always required, even for low-dim data
                    }
                }).then(
                    function(response) {
                        if (response.result.artifacts.files.length > 0) {
                            service.composeSummaryResults(response.result.artifacts.files, response.executionId, phase)
                                .then(
                                    function(result) { resolve({result: result}); },
                                    function(msg) { reject(msg.statusText); }
                                );
                        } else {
                            resolve({result: {}});
                        }
                    },
                    function(response) { reject(response); }
                );
            });
        };

        service.composeSummaryResults = function(files, executionId, phase) {
            // FIXME: errors from downloadJsonFile do not lead to a reject
            return $q(function(resolve, reject) {
                var retObj = {summary: [], allSamples: 0, numberOfRows: 0},
                    fileExt = {fetch: ['.png', 'json'], preprocess:['all.png', 'all.json']},

                // find matched items in an array by key
                    _find = function composeSummaryResults_find (key, array) {
                        // The variable results needs var in this case (without 'var' a global variable is created)
                        var results = [];
                        for (var i = 0; i < array.length; i++) {
                            if (array[i].search(key) > -1) {
                                results.push(array[i]);
                            }
                        }
                        return results;
                    },

                // process each item
                    _processItem  = function composeSummaryResults_processItem(img, json) {
                        return $q(function(resolve) {
                            service.downloadJsonFile(executionId, json).then(
                                function (d) {
                                    retObj.subsets = d.data.length;
                                    d.data.forEach(function(subset) {
                                        retObj.allSamples += subset.numberOfSamples;
                                        retObj.numberOfRows = subset.totalNumberOfValuesIncludingMissing /
                                            subset.numberOfSamples;
                                    });
                                    resolve({img: service.urlForFile(executionId, img), json:d});
                                },
                                function (err) { reject(err); }
                            );
                        });
                    };

                // first identify image and json files
                var _images = _find(fileExt[phase][0], files),
                    _jsons = _find(fileExt[phase][1], files);

                // load each json file contents
                for (var i = 0; i < _images.length; i++){
                    retObj.summary.push(_processItem(_images[i], _jsons[i]));
                }

                $.when.apply($, retObj.summary).then(function () {
                    resolve(retObj); // when all contents has been loaded
                });
            });
        };

        service.preprocess = function(args) {
            return $q(function (resolve, reject) {
                service.startScriptExecution({
                    taskType: 'preprocess',
                    arguments: args
                }).then(
                    resolve,
                    function(response) { reject(response); }
                );
            });
        };

        return service;
    }]);
