//# sourceURL=ppmidemo.js

'use strict';

window.smartRApp.controller('PPMIDemoController',
    ['$scope', 'smartRUtils', 'commonWorkflowService', function($scope, smartRUtils, commonWorkflowService) {

        commonWorkflowService.initializeWorkflow('ppmidemo', $scope);

        $scope.variantDB = {
            data: [],
            regions: '',
            genes: 'TP63',
            server: "http://bio3.uni.lu/accessDB/accessDB",
        };

        $scope.pdmap = {
            user: "test",
            password: "test.123",
            server: "http://pg-sandbox.uni.lu/minerva/galaxy.xhtml",
            model: 'pdmap_dec15',
        };

        var _emulateGridAnalysis = function() {
            var dfd = $.Deferred();

            runAllQueries(function() {
                gridstore = new Ext.data.JsonStore({
                    url: pageInfo.basePath+'/chart/analysisGrid'
                });

                gridstore.on('load', storeLoaded);

                var myparams = Ext.urlEncode({
                    concept_key: "",
                    result_instance_id1: GLOBAL.CurrentSubsetIDs[1],
                    result_instance_id2: GLOBAL.CurrentSubsetIDs[2]
                });

                gridstore.load({
                    params: myparams,
                    callback: function () {
                        dfd.resolve();
                    }
                });
            });

            return dfd.promise();
        }

        var getTMIDs = function() {
            var dfd = $.Deferred();
            runAllQueries(function() {
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
            return dfd.promise();
        };

        var getVariantDBIDs = function(tmIDs) {
            var path = '/individuals/POST/';
            var filter_command = 'getfields!eq!id,comments&comments!isa!' + tmIDs.map(function(d) { return d.id; }).join(',');
            return $.ajax({
                url: pageInfo.basePath + '/SmartR/variantDB',
                type: 'POST',
                data: {
                    server: $scope.variantDB.server,
                    path: path,
                    filter_command: filter_command,
                }
            }).then(
                function(res) {
                    return JSON.parse(res).values.map(function(d) { return d[0]; });
                },
                function() { alert('Connection refused: ' + $scope.variantDB.server); }
            );
        };

        var getVariantDBRequestsForGenes = function(ids, genes) {
            var path = '/variant_individuals/POST/';
            var filter_command = 'splitcommand!eq!t&getfields!eq!start,reference,alleleseq,variant_genotypes&shown_individuals!eq!' + ids.join(',') +
                '&variant_genotypes!ov!' + ids.join(',') + '&gene!eq!' + genes.join(',');
            return $.ajax({
                url: pageInfo.basePath + '/SmartR/variantDB',
                type: 'POST',
                data: {
                    server: $scope.variantDB.server,
                    path: path,
                    filter_command: filter_command,
                }
            }).then(
                function(res) { return JSON.parse(res); },
                function() { alert('Connection refused: ' + $scope.variantDB.server); }
            );
        }

        var getVariantDBData = function(request) {
            var path = '/variant_individuals/POST/';
            var filter_command = request;
            return $.ajax({
                url: pageInfo.basePath + '/SmartR/variantDB',
                type: 'POST',
                data: {
                    server: $scope.variantDB.server,
                    path: path,
                    filter_command: filter_command,
                }
            }).then(
                function(res) { return res; },
                function() { alert('Connection refused: ' + $scope.variantDB.server); }
            );
        };

        $scope.fetchVariantDB = function() {
            $scope.variantDB.data = [];
            var genes = $scope.variantDB.genes.split(',').map(function(d) { return d.trim(); });
            getTMIDs().then(function(tmIDs) {
                getVariantDBIDs(tmIDs).then(function(variantDBIDs) {
                    if ($scope.variantDB.genes) {
                        getVariantDBRequestsForGenes(variantDBIDs, genes).then(function(requests) {
                            requests.forEach(function(request) {
                                getVariantDBData(request).then(function(data) {
                                    var variantData = JSON.parse(data);
                                    var indices = {};
                                    indices['pos'] = variantData.fields.indexOf('start');
                                    indices['ref'] = variantData.fields.indexOf('reference');
                                    indices['alt'] = variantData.fields.indexOf('alleleseq');
                                    indices['chr'] = variantData.fields.indexOf('chrom');
                                    indices['ids'] = [];
                                    variantDBIDs.forEach(function(variantDBID) {
                                        var idx = variantData.fields.indexOf(variantDBID);
                                        if (idx !== -1) {
                                            indices['ids'].push(idx);
                                        }
                                    });

                                    variantData.values.forEach(function(d) {
                                        var pos = d[indices['pos']]
                                        var ref = d[indices['ref']]
                                        var alt = d[indices['alt']]
                                        var chr = d[indices['chr']]
                                        var variants = indices['ids'].map(function(idIndex) {
                                            return d[idIndex];
                                        });
                                        var frq = variants.filter(function(d) { return d.indexOf('1') !== -1; }) / variants.length;
                                        if (isNaN(frq)) { frq = 0; }

                                        $scope.variantDB.data.push({
                                            pos: pos,
                                            ref: ref,
                                            alt: alt,
                                            chr: chr,
                                            frq: frq,
                                            subset: 1, // FIXME: make dynamic
                                        });
                                    });
                                });
                            });
                        });
                    } else {

                    }
                })
            });
        };

        $scope.createPDMapLayout = function() {
            var expression_value = '#TYPE=GENETIC_VARIANT\n';
            expression_value += '#GENOME_TYPE=UCSC\n';
            expression_value += '#GENOME_VERSION=hg38\n';
            expression_value += 'position\toriginal_dna\talternative_dna\tname\tcontig\tallel_frequency\n';

            $scope.variantDB.data.forEach(function(d) {
                expression_value += d.pos + '\t' + d.ref + '\t' + d.alt + '\t' + 'FIXME' + '\t' + d.chr + '\t' + d.frq + '\n';
            });

            return $.ajax({
                url: pageInfo.basePath + '/SmartR/pdMap',
                type: 'POST',
                data: {
                    url: $scope.pdmap.server,
                    identifier: Math.random() * Math.pow(10,17),
                    login: $scope.pdmap.user,
                    password: $scope.pdmap.password,
                    model: $scope.pdmap.model,
                    expression_value: expression_value,
                }
            }).then(
                function(res) { console.log(res); },
                function(res) { console.error(res); }
            );
        };

    }]);

