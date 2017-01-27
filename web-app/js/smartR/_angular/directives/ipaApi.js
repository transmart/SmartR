//# sourceURL=ipaApi.js

'use strict';

angular.module('ipaApi',['ngTable'])
    .directive('ipaSimpleTable', [
        '$rootScope',
        function($rootScope) {
        return {
            restrict: 'E',
            scope: {
                simpleTable: '=',
                id: '=',
                tooltip: '=',
                numericPrecision: '=',
                paginationNumRows: '=',
            },
            css: [$rootScope.smartRPath + '/css/ng-table.css'],

            controller: ['$scope', 'NgTableParams', function($scope,NgTableParams) {
                function simpleTableToNgTable(simpleTable,precision) {
                    var ngtbl = {
                        header: [],
                        rows: [],
                    };
                    if (simpleTable == null)
                        return ngtbl;

                    // create column hashes as needed by ng-table
                    for (var c = 0; c < simpleTable.header.length; c++) {
                        var colname = simpleTable.header[c];
                        // sanitize field id FIXME is there a function to do that? Are we missing some letters?
                        var colid = colname.replace(/ /g,'')
                        colid = colid.replace(/#/g,'__hash__')
                        colid = colid.replace(/-/g,'__dash__')
                        colid = colid.replace(/\(/g,'__opnbr__')
                        colid = colid.replace(/\)/g,'__clsbr__')
                        ngtbl.header.push({ "title": colname, "field": colid, "visible": true, sortable: colid, });
                    }

                    // create one hash per row as needed by ng-table
                    for (var r = 0; r < simpleTable.rows.length; r++) {
                        var newrow = {};
                        for (var c = 0; c < ngtbl.header.length; c++) {
                            if (simpleTable.columnTypes != null && simpleTable.columnTypes[c] == "NUMERIC") {
                                // ensure that strings are converted to numbers, so sorting works
                                var thenum;
                                if (!simpleTable.rows[r][c] || /^\s*$/.test(simpleTable.rows[r][c])) {
                                    // empty string or only white space, set to NaN as it will be zero if we convert it using Number
                                    thenum = NaN;
                                } else {
                                    thenum = Number(simpleTable.rows[r][c]);
                                    if (isNaN(thenum) && typeof simpleTable.rows[r][c] == "string" && simpleTable.rows[r][c].toLowerCase() != "nan") {
                                        // warn if this was a string but not 'NaN'
                                        console.log("col "+simpleTable.header[c]+" ("+c+") NaN: "+simpleTable.rows[r][c]);
                                    }
                                }
                                if (precision != null) {
                                    // if a numeric precision was given, apply this
                                    thenum = Number(thenum.toPrecision(precision));
                                }
                                newrow[ngtbl.header[c].field] = thenum;
                            } else {
                                newrow[ngtbl.header[c].field] = simpleTable.rows[r][c];
                            }
                        }
                        ngtbl.rows.push(newrow);
                    }

                    return ngtbl;
                }

                function updateTable() {
                    var ngtbl = simpleTableToNgTable($scope.simpleTable,$scope.numericPrecision);
                    $scope.header = ngtbl.header;

                    var numRows = 20;
                    if ($scope.paginationNumRows != null) {
                        numRows = $scope.paginationNumRows;
                    }

                    // for table configuration see https://github.com/esvit/ng-table/wiki/Configuring-your-table-with-ngTableParams
                    $scope.tableParams = new NgTableParams(
                        {
                            count: numRows,
                        },
                        {
                            dataset: ngtbl.rows,
                            counts: [],
                        }
                    )
                }

                $scope.$watch('simpleTable', function(newValue) {
                    updateTable();
                }, true);

            }],

            template:
            '<table ng-table-dynamic="tableParams with header" class="table table-condensed table-bordered table-striped">' +
            '    <tr ng-repeat="row in $data">' +
            '    <td ng-repeat="col in $columns">{{row[col.field]}}</td>' +
            '    </tr>' +
            '</table>'
        }
    }])

    .directive('ipaApi', [
        '$rootScope',
        '$http',
        function($rootScope,$http) {
        return {
            restrict: 'E',
            scope: {
                differentiallyExpressed: '=',
                degParams: '=',
                credentials: '=',
            },

            controller: ['$scope', '$interval', '$filter', function($scope, $interval, $filter) {
                // intitialisation
                clearIpaState();

                function clearIpaState() {
                    $scope.ipaResults = null;
                    // this will trigger stopWaitingForIpaAnalysisId, which also resets launchIpaMsg
                    $scope.ipaState = {
                        analysisName: null,
                    };
                }

                // when differentiallyExpressed changed the IPA results are not valid anymore
                $scope.$watch('differentiallyExpressed', function(newValue) {
                    clearIpaState();
                }, true);

                function simpleTableToLists(simpleTable) {
                    var lists = {};
                    if (simpleTable == null)
                        return lists;

                    // initialise one list per column
                    for (var c = 0; c < simpleTable.header.length; c++) {
                        lists[simpleTable.header[c]] = [];
                    }

                    // create one hash per row as needed by ng-table
                    for (var r = 0; r < simpleTable.rows.length; r++) {
                        for (var c = 0; c < simpleTable.header.length; c++) {
                            lists[simpleTable.header[c]].push(simpleTable.rows[r][c]);
                        }
                    }

                    return lists;
                }

                function handleIpaFailure(response) {
                    console.log('Error at IPA backend, '+ response.status +' '+ response.statusText +': '+ response.data);
                    var longmsg;
                    if (response.data.length > 120) {
                        longmsg = response.data.substring(0,120)+' ...';
                    } else {
                        longmsg = response.data;
                    }
                    $scope.ipaState.launchIpaMsg = 'Error at IPA backend, '+ response.status +' '+ response.statusText +': '+ longmsg;
                    $scope.ipaState.launchIpaMsgBlink = false;
                };

                function lauchIpa() {
                    // convert differentiallyExpressed to the format needed by dataAnalysis
                    var tbllist = simpleTableToLists($scope.differentiallyExpressed);

                    // set the main data that will influence the results
                    var data = {
                        geneIdType: 'humanegsym', // Entrez Gene symbols
                        geneId: tbllist.symbol,
                        expValueType: 'logratio',
                        expValue: tbllist.logFC,
                        expValueType2: 'pvalue',
                    };
                    if ($scope.degParams.significanceMeasure == "pval") {
                        data.expValue2 = tbllist.pval;
                    } else if ($scope.degParams.significanceMeasure == "adjpval") {
                        data.expValue2 = tbllist.adjPval;
                    } else {
                        throw "unhandled significanceMeasure " + $scope.degParams.significanceMeasure
                    }

                    // generate a reproducible object (with defined order of keys)
                    var hashme = [];
                    Object.keys(data)
                        .sort()
                        .forEach(function(v, i) {
                            hashme.push(v);
                            hashme.push(data[v]);
                        });

                    // generate a reproducible hash of the data
                    var sha256hashed = Sha256.hash(JSON.stringify(hashme));

                    data.analysisName = sha256hashed;
                    data.applicationName = 'IPAtranSMART';
                    data.username = $scope.credentials.username
                    data.password = $scope.credentials.password

                    var url = pageInfo.basePath + '/IpaConnector/dataAnalysis';
                    var request = $http({
                        url: url,
                        method: 'POST',
                        config: {
                            timeout: 10000
                        },
                        data: data,
                    });

                    $scope.ipaState.launchIpaMsg = 'Launching Ingenuity Pathway Analysis'
                    $scope.ipaState.launchIpaMsgBlink = true;

                    request.then(
                        function(response) {
                            $scope.ipaState.launchIpaMsg = '';
                            $scope.ipaState.launchIpaMsgBlink = false;
                            // set the name *after* the call succeeded, because this will trigger polling the server
                            $scope.ipaState.analysisName = sha256hashed;
                        },
                        handleIpaFailure);
                }
                $scope.lauchIpa = lauchIpa;

                $scope.$watch('ipaState.analysisName', function(newValue) {
                    stopWaitingForIpaAnalysisId(); // this will also reset launchIpaMsg
                    if ($scope.ipaState.analysisName != null) {
                        $scope.ipaState.launchIpaMsg = 'Waiting for Ingenuity Pathway Analysis';
                        $scope.ipaState.launchIpaMsgBlink = true;
                        // start a timer to periodically check for analysis id
                        waitForIpaAnalysisId()
                        // try to get ID immediately (if it succeeds it will stop above timer)
                        analysisId($scope.ipaState.analysisName);
                    }
                }, true);

                var waitForIpaAnalysisIdTimer;
                function waitForIpaAnalysisId() {
                    // Don't start a new timer if it is running already
                    if ( angular.isDefined(waitForIpaAnalysisIdTimer) ) return;

                    waitForIpaAnalysisIdTimer = $interval(function() {
                        var curTime = $filter('date')(new Date(), 'HH:mm:ss');
                        $scope.ipaState.launchIpaMsg = 'Waiting for Ingenuity Pathway Analysis (last poll '+curTime+')';
                        $scope.ipaState.launchIpaMsgBlink = true;
                        analysisId($scope.ipaState.analysisName);
                    }, 20000);
                };

                function analysisId(analysisName) {
                    var data = {
                        analysisName: analysisName,
                        username: $scope.credentials.username,
                        password: $scope.credentials.password
                    };
                    var url = pageInfo.basePath + '/IpaConnector/analysisId';
                    var request = $http({
                        url: url,
                        method: 'POST',
                        config: {
                            timeout: 10000
                        },
                        data: data,
                    });

                    request.then(
                        function(response) {
                            if (response.data.analysisId != null)
                                $scope.ipaState.analysisId = response.data.analysisId;
                        },
                        handleIpaFailure);
                }

                $scope.$watch('ipaState.analysisId', function(newValue) {
                    stopWaitingForIpaAnalysisId();
                    if ($scope.ipaState.analysisId != null) {
                        exportIngenuityResults($scope.ipaState.analysisId);
                    }
                }, true);

                // cancel querying IPA server for results
                function stopWaitingForIpaAnalysisId() {
                    if (angular.isDefined(waitForIpaAnalysisIdTimer)) {
                        $interval.cancel(waitForIpaAnalysisIdTimer);
                        waitForIpaAnalysisIdTimer = undefined;
                    }
                    $scope.ipaState.launchIpaMsg = ''
                    $scope.ipaState.launchIpaMsgBlink = false;
                }

                function exportIngenuityResults(analysisId) {
                    var data = {
                        analysisId: analysisId,
                        username: $scope.credentials.username,
                        password: $scope.credentials.password
                    };
                    var url = pageInfo.basePath + '/IpaConnector/exportIngenuityResults';
                    var request = $http({
                        url: url,
                        method: 'POST',
                        config: {
                            timeout: 10000
                        },
                        data: data,
                    });

                    $scope.ipaState.launchIpaMsg = 'Exporting Ingenuity Pathway Analysis';
                    $scope.ipaState.launchIpaMsgBlink = true;
                    request.then(
                        function(response) {
                            $scope.ipaResults = response.data;
                        },
                        handleIpaFailure);
                }

                $scope.$on('$destroy', function() {
                    // make sure the interval is destroyed
                    stopWaitingForIpaAnalysisId();
                });

                $scope.$watch('ipaResults', function(newValue) {
                    $scope.ipaState.launchIpaMsg = '';
                    $scope.ipaState.launchIpaMsgBlink = false;
                    // reset tables
                    $scope.canonicalPathways = null;
                    $scope.upstreamRegulators = null;
                    $scope.toxFunctions = null;
                    $scope.networks = null;
                    $scope.diseasesAndBioFunctions = null;
                    $scope.analysisReadyMolecules = null;
                    if ($scope.ipaResults == null) {
                        // ensure that the Limma tab is shown when there are no IPA results;
                        $scope.ipatab = 1;
                    } else {
                        if ('tables' in $scope.ipaResults) {
                            if ('Canonical Pathways' in $scope.ipaResults.tables)
                                $scope.canonicalPathways = $scope.ipaResults.tables['Canonical Pathways'];
                            if ('Upstream Regulators' in $scope.ipaResults.tables)
                                $scope.upstreamRegulators = $scope.ipaResults.tables['Upstream Regulators'];
                            if ('Diseases and Bio Functions' in $scope.ipaResults.tables)
                                $scope.diseasesAndBioFunctions = $scope.ipaResults.tables['Diseases and Bio Functions'];
                            if ('Tox Functions' in $scope.ipaResults.tables)
                                $scope.toxFunctions = $scope.ipaResults.tables['Tox Functions'];
                            if ('Networks' in $scope.ipaResults.tables)
                                $scope.networks = $scope.ipaResults.tables['Networks'];
                            if ('Analysis Ready Molecules' in $scope.ipaResults.tables)
                                $scope.analysisReadyMolecules = $scope.ipaResults.tables['Analysis Ready Molecules'];
                        }
                    }
                }, true);
            }],

            // FIXME is there no other way to not depend on the $rootScope and smartR specifics
            templateUrl: $rootScope.smartRPath +  '/js/smartR/_angular/templates/ipaApi.html'
        }
        }]);
