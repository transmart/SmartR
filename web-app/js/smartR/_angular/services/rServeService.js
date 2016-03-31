//# sourceURL=rServeService.js

window.smartRApp.factory('rServeService', ['smartRUtils', '$q', '$http', function(smartRUtils, $q, $http) {

    var service = {};

    var NOOP_ABORT = function() {};
    var TIMEOUT = 10000 /* 10 s */;
    var CHECK_DELAY = 1000;
    var SESSION_TOUCH_DELAY = 9 * 60 * 1000; /* 9 min; session timeout is 10 */

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
        return $http({
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
        })
            .then(function(response) {
                state.sessionId = response.data.sessionId;
                rServeService_scheduleTouch();
            }, transformAjaxFailure);

    };

    service.touch = function(sessionId) {
        if (sessionId != state.sessionId) {
            return;
        }

        return $http({
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
        })
            .finally(function() {
                rServeService_scheduleTouch(); // schedule another
            })
    };

    function rServeService_scheduleTouch() {
        var sessionId = state.sessionId;
        window.clearTimeout(state.touchTimeout);
        state.touchTimeout = window.setTimeout(function() {
            service.touch(sessionId);
        }, SESSION_TOUCH_DELAY);
    }

    service.destroySession = function(sessionId) {
        sessionId = sessionId || state.sessionId;

        if (!sessionId) {
            return;
        }

        return $http({
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
        })
            .catch(transformAjaxFailure)
            .finally(function() {
                if (state.sessionId == sessionId) {
                    service.abandonCurrentSession();
                }
            });
    };

    service.abandonCurrentSession = function() {
        window.clearTimeout(state.touchTimeout);
        state.sessionId = null;
    };

    service.destroyAndStartSession = function(workflowName) {

        // delete session before creating a new one
        $q.when(service.destroySession())
            .then(function() {
                // start a new session
                service.startSession(workflowName);
            });

    };

    /*
     * taskData = {
     *     arguments: { ... },
     *     taskType: 'fetchData' or name of R script minus .R,
     *     phase: 'fetch' | 'preprocess' | 'run',
     * }
     *
     * If it succeeds, the result will be the data returned by the server.
     * If it fails, it will return an object with at least these fields:
     * {
     *   status: 0 | <http error code>,
     *   statusText: <a description of the error>,
     * }
     */
    service.startScriptExecution = function(taskDataOrig) {

        var taskData = $.extend({}, taskDataOrig); // clone the thing
        state.currentRequestAbort();

        var runRequest = $http({
            url: pageInfo.basePath + '/ScriptExecution/run',
            method: 'POST',
            timeout: TIMEOUT,
            responseType: 'json',
            data: JSON.stringify({
                sessionId: state.sessionId,
                arguments: taskData.arguments,
                taskType: taskData.taskType,
                workflow: workflow
            })
        });

        _setCancellationForAjaxCall(runRequest);

        /* schedule checks */
        var promise = runRequest.then(
            function(response) {
                taskData.executionId = response.data.executionId;
                return _checkStatus(taskData.executionId, CHECK_DELAY);
            },
            transformAjaxFailure
        );

        promise.cancel = function timeoutRequest_cancel() {
            // calling this method should by itself resolve the promise
            state.currentRequestAbort();
        };

        // no touching necessary when a task is running
        window.clearTimeout(state.touchTimeout);
        promise.finally(rServeService_scheduleTouch.bind(this));

        return promise;
    };

    function _setCancellationForAjaxCall(ajax) {
        /* request in-flight; aborting is cancelling this request */
        state.currentRequestAbort = function() { ajax.abort(); };
        /* once the request finishes, there's nothing to abort.
         * this needs to be the 1st callback, so that later callbacks can
         * override this */
        ajax.finally(function() {
            state.currentRequestAbort = NOOP_ABORT;
        });
    }

    /* set a function to be executed after a period of time and
     * simultaneously define the request cancellation function
     * as the cancellation of this timeout.
     * Return the result as a promise. */
    function _setStatusRequestTimeout(funcToCall, delay /*, ... */) {
        var defer = $.Deferred();
        var promise = defer.promise();
        var restOfArguments = Array.prototype.slice.call(arguments)
            .splice(2); // arguments after delay

        function _setStatusRequestTimeout_wrappedFunc() {
            // we cannot abort by calling clearTimeout() anymore
            // funcToCall will probably set its own abort method
            state.currentRequestAbort = NOOP_ABORT;

            // funcToCall should itself return a promise or a final result
            // but let's cover the case where it throws
            var result;
            try {
                result = funcToCall.apply(undefined, restOfArguments); // promise
                result.done(function(x) { defer.resolve(x); });
                result.fail(function(x) { defer.reject(x); });
            } catch (e) {
                defer.fail(e);
            }
        }

        var timeout = window.setTimeout(
            _setStatusRequestTimeout_wrappedFunc, delay);

        state.currentRequestAbort = function() {
            clearTimeout(timeout);
            if (defer.state() == 'pending') {
                defer.reject({
                    status: 0,
                    statusText: 'abort'
                });
            }
            state.currentRequestAbort = NOOP_ABORT;
        };

        return promise;
    }

    /* aux function of _startScriptExecution. Needs to follow its contract
     * with respect to the fail and success result of the promise */
    function _checkStatus(executionId, delay) {
        var ajax = $.ajax({
            type: 'GET',
            url : pageInfo.basePath + '/ScriptExecution/status',
            data: {
                sessionId  : state.sessionId,
                executionId: executionId
            }
        });

        state.currentRequestAbort = function() { ajax.abort(); };
        ajax.always(function() { state.currentRequestAbort = NOOP_ABORT; });

        return ajax.then(function (d) {
            if (d.state === 'FINISHED') {
                d.executionId = executionId;
                return d;
            } else if (d.state === 'FAILED') {
                return $.Deferred().reject({
                    status: 0,
                    statusText: d.result.exception
                }).promise();
            } else {
                // else still pending
                return _setStatusRequestTimeout(
                    _checkStatus, delay, executionId, delay);
            }
        }, transformAjaxFailure);
    }

    function transformAjaxFailure(response) {
        var ret = {
            status: response.status,
            statusText: response.statusText
        };
        if (response.data !== undefined) {
            ret.response = response.data;
        }

        return ret;
    }

    service.downloadJsonFile = function(executionId, filename) {
        // Simple GET request example:
        return $http({
            method: 'GET',
            url: this.urlForFile(executionId, filename)
        });
    };


    service.urlForFile = function(executionId, filename) {
        return pageInfo.basePath +
            '/ScriptExecution/downloadFile?sessionId=' +
            state.sessionId +
            '&executionId=' +
            executionId +
            '&filename=' +
            filename;
    };

    service.loadDataIntoSession = function(conceptKeys, dataConstraints) {
        return $q( function(resolve, reject) {
            smartRUtils.getSubsetIds().then(
                function(subsets) {
                    var _arg = {
                        conceptKeys: conceptKeys,
                        resultInstanceIds: subsets,
                        projection:'log_intensity'
                    };

                    if (typeof dataConstraints !== 'undefined') {
                        _arg.dataConstraints = dataConstraints;
                    }

                    service.startScriptExecution({
                        taskType: 'fetchData',
                        arguments: _arg
                    }).then(
                        function(ret) { resolve('Task complete! State: ' + ret.state); },
                        function(ret) { reject(ret.response); }
                    );
                },
                function() {
                    reject('Could not create subsets. Did you select a cohort?');
                }
            );
        });
    };

    service.executeSummaryStats = function(phase) {
        return $q( function(resolve, reject) {
            service.startScriptExecution({
                taskType: 'summary',
                arguments: {
                    phase: phase,
                    projection: 'log_intensity' // always required, even for low-dim data
                }
            }).then(
                function(ret) {
                    var _result = {};
                    if (ret.result.artifacts.files.length > 0) {
                        service.composeSummaryResults(ret.result.artifacts.files, ret.executionId, phase)
                            .then(function (result) {
                                _result = result;
                                resolve({result : _result, msg:'Task complete! State: ' + ret.state});
                            });
                    } else {
                        resolve({result : _result, msg:'Task complete! State: ' + ret.state});
                    }
                },
                function(ret) { reject(ret.response); }
            );
        });
    };

    service.composeSummaryResults = function(files, executionId, phase) {
        return $q(function (resolve, reject) {
            var retObj = {summary : [], allSamples : 0},
                fileExt = {fetch : ['.png', 'json'], preprocess :['all.png', 'all.json']},

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
                    return $q(function (resolve, reject) {
                        service.downloadJsonFile(executionId, json).then(
                            function (d) {
                                retObj.subsets = d.data.length;
                                d.data.forEach(function (subset) {
                                    retObj.allSamples += subset.numberOfSamples;
                                });
                                resolve({img: service.urlForFile(executionId, img), json:d})
                            },
                            function (err) {reject(err);}
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
                function(ret) { resolve('Task complete! State: ' + ret.state); },
                function(ret) { reject(ret.response); }
            );
        });
    };

    return service;
}]);
