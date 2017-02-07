//# sourceURL=variantmap.js

'use strict';

window.smartRApp.controller('VariantMapController',
    ['$scope', 'smartRUtils', 'commonWorkflowService', function($scope, smartRUtils, commonWorkflowService) {

        commonWorkflowService.initializeWorkflow('variantmap', $scope);

        $scope.fetch = {
            disabled: false,
            running: false,
            loaded: false,
            conceptBoxes: {
                highDimensional: {concepts: [], valid: true},
                numeric: {concepts: [], valid: true},
                categoric: {concepts: [], valid: true}
            },
            selectedBiomarkers: [],
        };

        $scope.runAnalysis = {
            disabled: true,
            running: false,
            scriptResults: {},
            params: {
                variants: [],
            },
        };


        $scope.variantDB = {
            server: "http://10.79.2.77/accessDB",
        };

        $scope.pdmap = {
            user: "test",
            password: "test.123",
            server: "http://pg-sandbox.uni.lu/minerva",
            servlet: "/galaxy.xhtml",
            model: 'pdmap_dec15',
        };

        $scope.messages = {
            error: "",
        };

        $scope.$watchGroup(['fetch.running', 'runAnalysis.running'],
            function(newValues) {
                var fetchRunning = newValues[0],
                    runAnalysisRunning = newValues[1];

                // clear old results
                if (fetchRunning) {
                    $scope.runAnalysis.scriptResults = {};
                }

                // disable tabs when certain criteria are not met
                $scope.fetch.disabled = runAnalysisRunning;
                $scope.runAnalysis.disabled = fetchRunning || !$scope.fetch.loaded;
            }
        );

        var getTMIDs = function() {
            var dfd = $.Deferred();
            runAllQueries(function() {
                $.ajax({
                    url: pageInfo.basePath + '/chart/clearGrid',
                    method: 'POST',
                    data: {
                        charttype: 'cleargrid',
                    }
                }).then(function() {
                    $.ajax({
                        url: pageInfo.basePath + '/chart/analysisGrid',
                        type: 'POST',
                        data: {
                            concept_key: '',
                            result_instance_id1: GLOBAL.CurrentSubsetIDs[1],
                            result_instance_id2: GLOBAL.CurrentSubsetIDs[2]
                        }
                    }).then(function(res) {
                        var ids = [];
                        JSON.parse(res).rows.map(function(d) {
                            ids.push({id: d.patient, subject: d.subject, subset: d.subset === 'subset1' ? 1 : 2});
                        });
                        dfd.resolve(ids);
                    });

                });
            });
            return dfd.promise();
        };

        var getVariantDBIDs = function(tmIDs) {
            var dfd = $.Deferred();
            $.ajax({
                url: pageInfo.basePath + '/SmartR/variantDB',
                type: 'POST',
                data: {
                    server: $scope.variantDB.server,
                    path: '/individuals/POST/',
                    filter_command: 'getfields!eq!id,comments&comments!eqa!' + tmIDs.map(function(d) { return d.id; }).map(function(d) { return d.toLowerCase(); }).join(','),
                },
                success: function(res) {
                    var data = JSON.parse(res);
                    var vIDIdx = data.fields.indexOf('id');
                    var tmIDIdx = data.fields.indexOf('comments');
                    var ids = [];
                    data.values.forEach(function(d) {
                        var hits = tmIDs.filter(function(e) { return e.id === d[tmIDIdx]; });
                        hits.forEach(function(hit) {
                            ids.push({
                                vID: d[vIDIdx],
                                tmID: hit.id,
                                subject: hit.subject,
                                subset: hit.subset
                            });
                        });
                    });
                    if (ids.length) {
                        dfd.resolve(ids);
                    } else {
                        dfd.reject('No matching Subject IDs found in VariantDB.');
                    }
                },
                failure: function() { dfd.reject('An error occured when trying to communicate with VariantDB.'); }
            });
            return dfd.promise();
        };

        var getVariantDBRequestsForGenes = function(ids) {
            var variantIDs = ids.map(function(d) { return d.vID; });
            var dfd = $.Deferred();
            var genes = $scope.fetch.selectedBiomarkers.map(function(d) { return d.name.toLowerCase(); });

            $.ajax({
                url: pageInfo.basePath + '/SmartR/variantDB',
                type: 'POST',
                data: {
                    server: $scope.variantDB.server,
                    path: '/variant_all/POST/',
                    filter_command: 'splitcommand!eq!t&getfields!eq!gene_at_position,start,reference,alleleseq,variant_genotypes,hgmd_rsid&shown_individuals!eq!' +
                    variantIDs.join(',') + '&c01,c11!ov!' + variantIDs.join(',') + '&gene!eq!' + genes.join(',')
                },
                success: function(res) { dfd.resolve(JSON.parse(res)); },
                failure: function() { dfd.reject('An error occured when trying to communicate with VariantDB.'); }

            });
            return dfd.promise();
        };

        var getVariantDBData = function(request) {
            var dfd = $.Deferred();
            $.ajax({
                url: pageInfo.basePath + '/SmartR/variantDB',
                type: 'POST',
                data: {
                    server: $scope.variantDB.server,
                    path: '/variant_all/POST/',
                    filter_command: request,
                },
                success: function(res) {
                    dfd.resolve(res);
                },
                failure: function() {
                    dfd.reject('An error occured when trying to communicate with VariantDB.');
                }
            });
            return dfd.promise();
        };

        var _prepareData = function(data, subset, subsetIDs) {
            var variantData = JSON.parse(data);
            var indices = {};
            indices['pos'] = variantData.fields.indexOf('start');
            indices['ref'] = variantData.fields.indexOf('reference');
            indices['alt'] = variantData.fields.indexOf('alleleseq');
            indices['chr'] = variantData.fields.indexOf('chrom');
            indices['gene'] = variantData.fields.indexOf('gene_at_position');
            indices['rsid'] = variantData.fields.indexOf('hgmd_rsid');
            indices['ids'] = [];
            subsetIDs.forEach(function(id) {
                var idx = variantData.fields.indexOf(id.vID);
                if (idx !== -1) {
                    indices['ids'].push(idx);
                }
            });

            variantData.values.forEach(function(d) {
                var pos = d[indices['pos']];
                var ref = d[indices['ref']];
                var alt = d[indices['alt']];
                var chr = d[indices['chr']];
                var gene = d[indices['gene']];
                var rsid = d[indices['rsid']]
                var variants = indices['ids'].map(function(idIndex) {
                    var variant = d[idIndex];
                    var matches = variant.match(/1/g);
                    var count = matches || 0;
                    return (count === 2 ? 1 : 0) + '|' + (count > 0 ? 1 : 0);
                });
                var frq = variants.filter(function(d) { return d.indexOf(1) !== -1; }).length / subsetIDs.length;
                if (! isNaN(frq) && frq > 0) {
                    variants.forEach(function(variant, idx) {
                        if (variant.indexOf(1) === -1) {
                            return;
                        }
                        $scope.runAnalysis.params.variants.push({
                            subject: subsetIDs[idx].subject,
                            gene: gene,
                            pos: pos,
                            ref: ref,
                            alt: alt,
                            chr: chr,
                            var: variant,
                            rsid: rsid,
                            frq: frq,
                            subset: subset,
                        });
                    });
                }
            });
        };

        var prepareData = function(data, ids) {
            var subsets = smartRUtils.unique(ids.map(function(d) { return d.subset; }));
            subsets.forEach(function(subset) {
                var subsetIDs = ids.filter(function(d) { return d.subset === subset; });
                _prepareData(data, subset, subsetIDs);
            });
        };

        var handleRequests = function(requests, ids, dfd) {
            if (! requests.length) {
                dfd.resolve();
            } else {
                $.when(getVariantDBData(requests[0])).then(function(data) {
                    prepareData(data, ids);
                    handleRequests(requests.slice(1), ids, dfd);
                }, function(err) {
                    $scope.messages.error += '<br/>' + err;
                });
            }
        };

        var cleanup = function() {
            $scope.messages.error = '';
            $scope.runAnalysis.params.variants = [];
        };

        $scope.fetchVariantDB = function() {
            var dfd = $.Deferred();
            cleanup();
            if (!$scope.fetch.selectedBiomarkers.length) {
                dfd.reject('You must select one or more genes to continue.');
            } else {
                getTMIDs().then(function(tmIDs) {
                    getVariantDBIDs(tmIDs).then(function(ids) {
                        getVariantDBRequestsForGenes(ids).then(function(requests) {
                            $scope.$apply();
                            handleRequests(requests, ids, dfd);
                        }, dfd.reject);
                    }, dfd.reject);
                });
            }
            return dfd.promise();
        };

        var _createPDMapLayout = function(subset, identifier) {
            var expression_value = '#TYPE=GENETIC_VARIANT\n';
            expression_value += '#GENOME_TYPE=UCSC\n';
            expression_value += '#GENOME_VERSION=hg38\n';
            expression_value += 'position\toriginal_dna\talternative_dna\tname\tcontig\tallel_frequency\n';

            $scope.variantDB.data.forEach(function(d) {
                if (d.subset === subset) {
                    expression_value += d.pos + '\t' + d.ref + '\t' + d.alt + '\t' + d.gene + '\t' + d.chr + '\t' + d.frq + '\n';
                }
            });

            return $.ajax({
                url: pageInfo.basePath + '/SmartR/pdMap',
                type: 'POST',
                data: {
                    url: $scope.pdmap.server + $scope.pdmap.servlet,
                    identifier: identifier,
                    login: $scope.pdmap.user,
                    password: $scope.pdmap.password,
                    model: $scope.pdmap.model,
                    expression_value: expression_value,
                }
            });
        };

        $scope.createPDMapLayout = function() {
            var cohorts = smartRUtils.countCohorts();
            var identifier1 = Math.random() * Math.pow(10,17);
            var identifier2 = Math.random() * Math.pow(10,17);
            _createPDMapLayout(1, identifier1).then(function() {
                if (cohorts === 2) {
                    _createPDMapLayout(2, identifier2).then(function() {
                        window.open($scope.pdmap.server);
                    });
                } else {
                    window.open($scope.pdmap.server);
                }
            });
        };

    }]);

