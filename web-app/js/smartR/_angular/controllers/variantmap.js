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
                method: 'pearson',
                transformation: 'raw'
            }
        };

        $scope.variantDB = {
            data: [],
            selectedGenes: '',
            server: "http://bio3.uni.lu/accessDB/accessDB",
            func_refgene: {
                exonic: false,
                intronic: false,
                upstream: false,
                downstream: false,
                splicing: false,
                intergenetic: false
            },
            exonicfunc_refgene: {
                frameshift_insertion: false,
                nonframeshift_insertion: false,
                frameshift_deletion: false,
                nonframeshift_deletion: false,
                frameshift_substitution: false,
                nonframeshift_substitution: false,
                synonymous_SNV: false,
                nonsynonymous_SNV: false,
            },
            misc: {
                globalMAF: 0.1
            },
            invalid: false,
            running: false
        };

        $scope.pdmap = {
            user: "test",
            password: "test.123",
            server: "http://pg-sandbox.uni.lu/minerva",
            servlet: "/galaxy.xhtml",
            model: 'pdmap_dec15',
            invalid: true,
        };

        $scope.intersection = {
            genes: []
        };

        $scope.messages = {
            error: "",
            loading: "",
            totalRequests: 0,
            finishedRequests: 0,
        };

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
                            ids.push({id: d.patient, subset: d.subset === 'subset1' ? 1 : 2});
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

        var getFilterString = function() {
            var filters1 = [];
            var filters2 = [];
            var filtersString = '';
            for (var key1 in $scope.variantDB.func_refgene) {
                if ($scope.variantDB.func_refgene.hasOwnProperty(key1) && $scope.variantDB.func_refgene[key1]) {
                    filters1.push(key1);
                }
            }
            for (var key2 in $scope.variantDB.exonicfunc_refgene) {
                if ($scope.variantDB.exonicfunc_refgene.hasOwnProperty(key2) && $scope.variantDB.exonicfunc_refgene[key2]) {
                    filters1.push(key2);
                }
            }

            filtersString += filters1.length ? '&func_refgene!ov!' + filters1.join(',') : filtersString;
            filtersString += filters2.length ? '&exonicfunc_refgene!ov!' + filters2.join(',') : filtersString;
            return filtersString;
        };

        var getVariantDBRequestsForGenes = function(variantDBIDs) {
            var ids = variantDBIDs.map(function(d) { return d.vID; });
            var dfd = $.Deferred();
            var genes = [];
            if ($scope.variantDB.selectedGenes.length) {
                genes = $scope.variantDB.selectedGenes.split(',').map(function(d) { return d.trim().toLowerCase(); })
                    .filter(function(d) { return $scope.intersection.genes.indexOf(d) !== -1; });
                if (! genes.length) {
                    dfd.reject("None of the specified genes could be found in VariantDB.");
                }
            } else {
                genes = $scope.intersection.genes;
            }

            $.ajax({
                url: pageInfo.basePath + '/SmartR/variantDB',
                type: 'POST',
                data: {
                    server: $scope.variantDB.server,
                    path: '/variant_all/POST/',
                    filter_command: 'splitcommand!eq!t&getfields!eq!gene_at_position,start,reference,alleleseq,variant_genotypes&shown_individuals!eq!' + ids.join(',') +
                        '&c01,c11!ov!' + ids.join(',') + '&gene!eq!' + genes.join(',') + getFilterString()
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

        var _prepareData = function(data, subset, ids) {
            var variantData = JSON.parse(data);
            var indices = {};
            indices['pos'] = variantData.fields.indexOf('start');
            indices['ref'] = variantData.fields.indexOf('reference');
            indices['alt'] = variantData.fields.indexOf('alleleseq');
            indices['chr'] = variantData.fields.indexOf('chrom');
            indices['gene'] = variantData.fields.indexOf('gene_at_position');
            indices['ids'] = [];
            ids.forEach(function(variantDBID) {
                var idx = variantData.fields.indexOf(variantDBID);
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
                var variants = indices['ids'].map(function(idIndex) {
                    return d[idIndex];
                });
                var frq = variants.filter(function(d) { return d.indexOf('1') !== -1; }).length / ids.length;
                if (! isNaN(frq) && frq !== 0 && frq >= $scope.variantDB.misc.cohortMAF) {
                    $scope.variantDB.data.push({
                        pos: pos,
                        ref: ref,
                        alt: alt,
                        chr: chr,
                        frq: frq,
                        gene: gene,
                        subset: subset,
                    });
                }
            });
        };

        var prepareData = function(data, variantDBIDs) {
            var subsets = smartRUtils.unique(variantDBIDs.map(function(d) { return d.subset; }));
            subsets.forEach(function(subset) {
                var ids = variantDBIDs.filter(function(d) { return d.subset === subset; })
                    .map(function(d) { return d.vID; });
                _prepareData(data, subset, ids);
            });
            checkPDMapSanity();
            $scope.$apply();
        };

        var handleError = function(err) {
            $scope.variantDB.running = false;
            $scope.messages.finishedRequests += 1;
            $scope.messages.error = err;
            $scope.$apply();
        };

        var handleRequests = function(requests, variantDBIDs) {
            if (! requests.length && $scope.messages.finishedRequests === $scope.messages.totalRequests) {
                $scope.variantDB.running = false;
                checkPDMapSanity();
                $scope.$apply();
                return;
            }
            if (! $scope.variantDB.running) { // if error
                return;
            }
            $.when(getVariantDBData(requests[0])).then(function(data) {
                prepareData(data, variantDBIDs);
                $scope.messages.finishedRequests += 1;
                $scope.$apply();
                handleRequests(requests.slice(1), variantDBIDs);
            }, handleError);
        };

        var cleanup = function() {
            $scope.variantDB.data = [];
            $scope.messages.error = '';
            $scope.messages.totalRequests = 0;
            $scope.messages.finishedRequests = 0;
            $scope.messages.loading = '';
        };

        $scope.fetchVariantDB = function() {
            cleanup();
            $scope.variantDB.running = true;
            getTMIDs().then(function(tmIDs) {
                $scope.variantDB.running = true;
                getVariantDBIDs(tmIDs).then(function(variantDBIDs) {
                    getVariantDBRequestsForGenes(variantDBIDs).then(function(requests) {
                        $scope.messages.totalRequests += requests.length;
                        $scope.$apply();
                        handleRequests(requests, variantDBIDs);
                    }, handleError);
                }, handleError);
            });
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

