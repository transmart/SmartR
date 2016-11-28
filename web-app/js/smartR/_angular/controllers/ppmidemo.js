//# sourceURL=ppmidemo.js

'use strict';

window.smartRApp.controller('PPMIDemoController',
    ['$scope', 'smartRUtils', 'commonWorkflowService', function($scope, smartRUtils, commonWorkflowService) {

        commonWorkflowService.initializeWorkflow('ppmidemo', $scope);

        $scope.variantDB = {
            data: [],
            regions: '',
            selectedGenes: '',
            server: "http://bio3.uni.lu/accessDB/accessDB",
            invalid: false,
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
        }

        $scope.messages = {
            error: "",
        };

        $scope.$watch('variantDB', function(vdb) {
            $scope.variantDB.invalid = !vdb.server || (!!vdb.regions && !!vdb.selectedGenes);
        }, true);

        var checkPDMapSanity = function() {
            $scope.pdmap.invalid = !$scope.pdmap.user || !$scope.pdmap.password || !$scope.pdmap.server || !$scope.variantDB.data.length;
        };
        $scope.$watch('pdmap', checkPDMapSanity, true);

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

        var setIntersectionGenes = function() {
            $.ajax({
                url: pageInfo.basePath + '/SmartR/variantDB',
                type: 'POST',
                async: false,
                data: {
                    server: $scope.variantDB.server,
                    path: '/validgenes/',
                    filter_command: ''
                },
                success: function(res) {
                    var variantDBGeneList = JSON.parse(res);
                    variantDBGeneList = variantDBGeneList.map(function(d) { return d.toLowerCase(); });

                    var PDMAP_GENES = ['LAG3','C19Orf12','C9Orf72','CHMP2B','DJ1','DNAJC13','FUS','GRN','PANK2','PRKN','SMPD1','SPG11','TARDBP','VCP','WDR45','TNK2','TNR','SLC39A8','DNM3','UCP2','SOD2','COL12A1','VAPB','ANKRD13A','MKS1','TOR1A','PTEN','ASNA1','PML','VPS53','SLC52A1','SLC5A9','EPPK1','FBXL17','RUNDC2A','KNDC1','CHCHD2','HNF4A','PTBP1','RAB39B','AAK1','ABCA1','ABCA5','ACMSD','AD7CNTP','ADAMTS16','ADAMTSL1','ADH1C','AGAP1','AIG1','AIMP1','AKAP13','AKT1','ANO5','AP1M1','ARRDC4','AS3MT','ASS1P9','ATF6','ATP12A','ATP13A2','ATP1A3','ATP6AP2','ATXN2','ATXN3','AXIN1','B3GALT2','BANF1P1','BANK1','BCL11A','BMP4','BRINP1','BST1','BTF3P2','BTNL2','C10orf32','C11orf71','C12orf55','C17orf51','C17orf97','C18orf1','C1QTNF6','C20orf85','C3orf20','C5orf64','C8orf83','C8orf86','C9','CADM2','CASC15','CAST','CCDC62','CCNY','CCT5P2','CD38','CD93','CDC14C','CDCA7L','CDCP2','CDH8','CDK17','CDK5RAP2','CEBPB','CENPC','CENPC1','CEP152','CERS6','CGRRF1','CHI3L2','CHIA','CHORDC1','CNKSR3','CNOT6','CNTNAP2','COL22A1','COL2A1','COL4A1','COLGALT2','COMT','COX6CP2','CREM','CRHBP','CRHR1','CRHR1-IT1','CSF1','CSMD1','CSTL1','CTIF','CTNNA3','CTSB','CUL2','CYB5RL','CYCSP38','CYCSP42','CYMP','CYP17A1','CYP2D6','CYP51A1P1','CYTL1','DAW1','DBC1','DCDC2C','DCTN1','DGCR5','DGCR9','DGKQ','DHFRP2','DIAPH3','DIS3L2P1','DKK2','DLG2','DNAJC1','DNAJC5','DNAJC6','DPPA5','DPY19L2P3','DSG3','DUSP27','EFCAB2','EFCAB4B','EIF3FP3','EIF4A1P6','EIF4G1','EPS8L3','ERMP1','ETV6','EVX2','EXD2','EXT1','FAM19A5','FAM32E','FAM47E','FAM49A','FAM69A','FBN1','FBXO3','FBXO7','FDFT1','FGF12','FLJ41941','FLJ45872','FOSL2','FOXF1','FRG1','FRMD4A','FRMPD4','FTL','GABBR2','GABRG3','GAK','GALNT16','GBA','GCH1','GFPT2','GIGYF2','GLT25D2','GPNMB','GPR156','GPRIN3','GRIN2A','GSTP1','GXYLT1','HACE1','HIP1R','HIVEP2','HLA-DR85','HLA-DRA','HLA-DRB5','HLA-S','HNRNPU','HSP90AA4P','HTRA2','HYI','IGLL1','IGSF5','IL1RL2','IL2RB','IMP5','INS-IGF2','IRF8','IRS1','ISCA1P2','ISM1','ITGA8','ITGAL','ITSN2','KANSL1','KBTBD11','KCNA10','KCNK9','KIAA0427','KIAA0947','KIAA1267','KIAA1279','KIAA1432','KIAA1715','KIF2A','KLF2','LAMP1','LAMP3','LCN1P2','LDHAP1','LDLRAD4','LEKR1','LINC00487','LINC00632','LMO2','LMOD1','LPHN2','LRRC37A','LRRK2','LTBP1','LYAR','MACROD2','MAF','MAP2K6','MAP3K1','MAP3K4','MAPK6PS2','MAPRE1P1','MAPRE2','MAPT','MBL2','MCCC1','MCHR2','MED13','MESTP4','MGAT5','MGC57346','MMRN1','MRPL50P1','MRPS35P3','MRPS36P2','MT1F','MT1G','NASPP1','NDUFS4','NEGR1','NELL1','NNMT','NOS1','NOS2A','NPM1P10','NPM1P18','NPY5R','NSF','NUB1','NUCKS1','NXT1','NYAP2','OCA2','OR7E28P','OR7E90P','OR9G3P','OR9G9','OTOG','PAICSP3','PAK1','PARK2','PARK7','PAX3','PAX7','PCDH15','PCDH18','PCDH20','PDZRN4','PEX26','PHACTR2','PHF5GP','PHKA2','PINK1','PITX2','PKDCC','PLA2G6','PLB1','PLCB4','PLEKHM1','PLG','PLXNC1','PM20D1','PMEPA1','POLG','POLR3B','POU6F2','PPEF1','PRDM13','PRDM15','PRDM2','PRKAG2','PRKRA','PRRG4','PRSS53','PTGS2','PTPRF','PTPRR','QSER1','RAB25','RAB29','RAD18','RAD51AP2','RAG1AP1','RAI1','RAP1A','RAPGEF5','RBFOX3','RBM26','RFX4','RGS2','RHOU','RIT2','RNF114','RNF5P1','RNU5B-6P','RPH3AL','RPL10AP7','RPL21P36','RPL23AP22','RPL23AP28','RPL26P19','RPL26P3','RPL34P11','RPL35AP12','RPL36AP23','RPL7AP57','RPL7P27','RPL9P21','RPS12P16','RPS12P4','RPS12P5','RPS12P8','RPS17P6','RPS19P6','RPS20P25','RPS23P3','RPS29P1','RPS2P34','RPS3AP46','RPS3P6','RPS3P7','RPS6P12','RREB1','RTN4R','RUFY3','SAMD12','SAMD4A','SCN2A','SCNN1A','SEMA5A','SEMA6D','SETD1A','SFRP1','SGCZ','SGOL1','SH3GL2','SIRT2','SLC17A6','SLC2A13','SLC41A1','SLC50A1','SLC6A3','SLC7A11','SNCA','SNORD50B','SNX10','SORCS1','SOX11','SOX3','SPATA8','SPHKAP','SPPL2C','SQRDL','SREBF1','SREBF2','ST13P2','STAP1','STH','STK39','SUN3','SYNJ1','SYT1','SYT11','TAF1','TAS1R2','TASP1','TBC1D22A','TBP','TCEAL2','TCF12','TFPI','TGFBR2','TH','THRB','TIMM17A','TLR4','TMEM108','TMEM128','TMEM132B','TMEM163','TMEM175','TMEM55A','TMEM72','TMPRSS6','TNF','TNFA','TNPO2','TPK1','TREM2','TRIM51FP','TRIQK','TRNAE39P','TRPM7','TTLL7','TYRL','TYW1B','UBBP4','UCHL1','UCHL5','ULK2','UNC13B','UNC13C','USP15','USP25','USP36','VENTXP7','VPS13C','VPS35','VPS41','VPS8','VWA5A','VWC2','WIPF3','WNT3','ZDHHC8','ZDHHC8P1','ZFP42','ZMAT4','ZNF280D','ZNF299P','ZNF385B','ZNF615','ZNF646','ZWINT'].map(function(d) { return d.toLowerCase(); });

                    var intersection = smartRUtils._intersectArrays(variantDBGeneList, PDMAP_GENES);
                    $scope.intersection.genes = intersection;
                }
            });
        };
        setIntersectionGenes();

        var getVariantDBIDs = function(tmIDs) {
            var dfd = $.Deferred();
            $.ajax({
                url: pageInfo.basePath + '/SmartR/variantDB',
                type: 'POST',
                data: {
                    server: $scope.variantDB.server,
                    path: '/individuals/POST/',
                    filter_command: 'getfields!eq!id,comments&comments!isa!' + tmIDs.map(function(d) { return d.id; }).join(','),
                },
                success: function(res) {
                    var ids = JSON.parse(res).values.map(function(d) { return d[0]; });
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
                        '&variant_genotypes!ov!' + ids.join(',') + '&gene!eq!' + genes.join(','),
                },
                success: function(res) { dfd.resolve(JSON.parse(res)); },
                failure: function() { dfd.reject('An error occured when trying to communicate with VariantDB.'); }

            });
            return dfd.promise();
        };

        var getVariantDBRequestsForRegions = function(ids) {
            var dfd = $.Deferred();
            var regions = $scope.variantDB.regions;
            regions = regions.split(',').map(function(d) { return d.trim(); });
            var chrs = regions.map(function(d) { return d.split(':')[0]; });
            var starts = regions.map(function(d) { return d.split(':')[1].split('-')[0]; });
            var stops = regions.map(function(d) { return d.split(':')[1].split('-')[1]; });
            $.ajax({
                url: pageInfo.basePath + '/SmartR/variantDB',
                type: 'POST',
                data: {
                    server: $scope.variantDB.server,
                    path: '/variant_all/POST/',
                    filter_command: 'splitcommand!eq!t&getfields!eq!start,reference,alleleseq,variant_genotypes&shown_individuals!eq!' +
                            ids.join(',') + '&variant_genotypes!ov!' + ids.join(',') + '&chrom!eq!' + chrs.join(',') + '&start!gt!' +
                            starts.join(',') + '&start!lt!' + stops.join(',')
                },
                success: function(res) { dfd.resolve(res); },
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
                success: function(res) { dfd.resolve(res); },
                failure: function() { dfd.reject('An error occured when trying to communicate with VariantDB.'); }
            });
            return dfd.promise();
        };

        var prepareData = function(data, subset, variantDBIDs) {
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
            $scope.$apply(function() {
                checkPDMapSanity();
            });
        };

        var handleError = function(err) {
            $scope.messages.error = err;
            $scope.$apply();
        };

        $scope.fetchVariantDB = function() {
            $scope.messages.error = '';
            $scope.variantDB.data = [];
            getTMIDs().then(function(tmIDs) {
                var subsets = smartRUtils.unique(tmIDs.map(function(d) { return d.subset; }));
                subsets.forEach(function(subset) {
                    var tmSubsetIDs = tmIDs.filter(function(d) { return d.subset === subset; });
                    getVariantDBIDs(tmSubsetIDs).then(function(variantDBIDs) {
                        if ($scope.variantDB.regions) {
                            getVariantDBRequestsForRegions(variantDBIDs).then(function(requests) {
                                requests.forEach(function(request) {
                                    getVariantDBData(request).then(function(data) {
                                        prepareData(data, subset, variantDBIDs);
                                    });
                                });
                            });
                        } else {
                            getVariantDBRequestsForGenes(variantDBIDs).then(function(requests) {
                                requests.forEach(function(request) {
                                    getVariantDBData(request).then(function(data) {
                                        prepareData(data, subset, variantDBIDs);
                                    }, handleError);
                                });
                            }, handleError);
                        }
                    }, handleError);
                });
            }, handleError);
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

