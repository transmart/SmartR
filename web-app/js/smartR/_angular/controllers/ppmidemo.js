//# sourceURL=ppmidemo.js

'use strict';

window.smartRApp.controller('PPMIDemoController',
    ['$scope', 'smartRUtils', 'commonWorkflowService', function($scope, smartRUtils, commonWorkflowService) {

        commonWorkflowService.initializeWorkflow('ppmidemo', $scope);

        $scope.variantDB = {
            data: [],
            regions: '',
            genes: 'PARK2',
            server: "http://bio3.uni.lu/accessDB/accessDB",
        };

        $scope.pdmap = {
            user: "test",
            password: "test.123",
            server: "http://pg-sandbox.uni.lu/minerva",
            servlet: "/galaxy.xhtml",
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
            var path = '/variant_all/POST/';
            var filter_command = 'splitcommand!eq!t&getfields!eq!gene_at_position,start,reference,alleleseq,variant_genotypes&shown_individuals!eq!' + ids.join(',') +
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
        };

        var getVariantDBData = function(request) {
            var path = '/variant_all/POST/';
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
                var subsets = smartRUtils.unique(tmIDs.map(function(d) { return d.subset; }));
                subsets.forEach(function(subset) {
                    var tmSubsetIDs = tmIDs.filter(function(d) { return d.subset === subset; });
                    getVariantDBIDs(tmSubsetIDs).then(function(variantDBIDs) {
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
                                        indices['gene'] = variantData.fields.indexOf('gene_at_position');
                                        indices['ids'] = [];
                                        variantDBIDs.forEach(function(variantDBID) {
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
                                            var frq = variants.filter(function(d) { return d.indexOf('1') !== -1; }).length / variants.length;
                                            if (isNaN(frq)) { frq = 0; }

                                            $scope.variantDB.data.push({
                                                pos: pos,
                                                ref: ref,
                                                alt: alt,
                                                chr: chr,
                                                frq: frq,
                                                gene: gene,
                                                subset: subset,
                                            });
                                        });
                                    });
                                });
                            });
                        } else {

                        }
                    })
                });
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
                        window.open($scope.pdmap.server + '?id=' + $scope.pdmap.model + '&layout=' + identifier1);
                        window.open($scope.pdmap.server + '?id=' + $scope.pdmap.model + '&layout=' + identifier2);
                    });
                } else {
                    window.open($scope.pdmap.server + '?id=' + $scope.pdmap.model + '&layout=' + identifier1);
                }
            });
        };

    }]);

