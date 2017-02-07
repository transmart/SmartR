//# sourceURL=d3VariantMap.js

'use strict';

window.smartRApp.directive('variantMap', [
    'smartRUtils',
    '$rootScope',
    function(smartRUtils, $rootScope) {

        return {
            restrict: 'E',
            scope: {
                data: '=',
                width: '@',
                height: '@'
            },
            templateUrl: $rootScope.smartRPath + '/js/smartR/_angular/templates/variantMap.html',
            link: function (scope, element) {
                var viz = element.children()[1];
                scope.$watch('data', function() {
                    $(viz).empty();
                    if (! $.isEmptyObject(scope.data)) {
                        smartRUtils.prepareWindowSize(scope.width, scope.height);
                        createVariantMap(scope, viz);
                    }
                });
            }
        };

        function createVariantMap(scope, viz) {
            var exprNum = smartRUtils.getElementWithoutEventListeners('sr-vm-expr-num');
            exprNum.addEventListener('input', function() {
                drawBoxes();
            });
            var cf = crossfilter(scope.data.data);
            var bySubject = cf.dimension(function(d) { return d.subject; });
            var byGene = cf.dimension(function(d) { return d.gene; });

            var subjects = smartRUtils.unique(getValuesForDimension(bySubject));
            var genes = smartRUtils.unique(getValuesForDimension(byGene));
            var clinicalFeatures = smartRUtils.unique(scope.data.data.reduce(function(prev, current) {
                return prev.concat(Object.keys(current));
            }, [])).filter(function(d) { return d.indexOf('fullName:') !== -1; });

            var BOX_HEIGHT = 30;
            var BOX_WIDTH = 10;
            var CLINICAL_BOX_HEIGHT = 10;

            var MAIN_PLOT_OFFSET = clinicalFeatures.length * CLINICAL_BOX_HEIGHT + 2;

            var height = genes.length * BOX_HEIGHT + MAIN_PLOT_OFFSET;
            var width = subjects.length * BOX_WIDTH;

            var margin = {
                top: 20,
                bottom: 200,
                right: 50,
                left: 100,
            };

            smartRUtils.prepareWindowSize(width + margin.left + margin.right, height + margin.top + margin.bottom);

            var svg = d3.select(viz).append('svg')
                .attr('width', width + margin.left + margin.right)
                .attr('height', height + margin.top + margin.bottom)
                .append('g')
                .attr('transform', 'translate(' + margin.left + ',' + margin.top + ')');

            function drawBoxes() {
                var boxData = [];
                subjects.forEach(function(subject) {
                    bySubject.filterExact(subject);
                    genes.forEach(function(gene) {
                        byGene.filterExact(gene);
                        var hits = bySubject.top(Infinity);
                        var zscore = hits.length ? hits[0].zscore : undefined;
                        var consequences = smartRUtils.unique(hits.reduce(function(prev, curr) {
                            if (curr.consequence) {
                                return prev.concat(curr.consequence.split(','));
                            } else {
                                return prev;
                            }
                        }, []));
                        var predictions = hits.map(function(d) { return d.prediction; });
                        var risk = (function() {
                            if (predictions.indexOf('probably damaging') !== -1) {
                                return 2;
                            } else if (predictions.indexOf('possibly damaging') !== -1) {
                                return 1;
                            } else if (predictions.indexOf('benign') !== -1) {
                                return 0;
                            } else {
                                return -1;
                            }
                        })();
                        var aachange = !!consequences.filter(function(consequence) { return consequence !== 'synonymous'; }).length;
                        boxData.push({
                            subject: subject,
                            gene: gene,
                            zscore: zscore,
                            consequences: consequences,
                            aachange: aachange,
                            variants: hits.length,
                            predictions: predictions,
                            risk: risk,
                        });
                        byGene.filterAll();
                    });
                    bySubject.filterAll();
                });
                // DATA JOIN
                var box = svg.selectAll('.sr-vm-box')
                    .data(boxData, function(d) { return d.subject + '-' + d.gene; });

                // ENTER g
                var boxEnter = box.enter()
                    .append('g')
                    .attr('class', 'sr-vm-box');

                // ENTER rect (main)
                boxEnter.append('rect')
                    .attr('class', 'sr-vm-main-box')
                    .attr('height', BOX_HEIGHT)
                    .attr('width', BOX_WIDTH);

                var mainBoxStroke = parseInt(d3.select('.sr-vm-main-box').style('stroke-width'));
                // ENTER rect (expr)
                boxEnter.append('rect')
                    .attr('class', 'sr-vm-expr-box')
                    .attr('height', BOX_HEIGHT - mainBoxStroke)
                    .attr('width', BOX_WIDTH - mainBoxStroke)
                    .attr('x', mainBoxStroke / 2)
                    .attr('y', mainBoxStroke / 2);

                // UPDATE g
                box.attr('transform', function(d) {
                    return 'translate(' + (subjects.indexOf(d.subject) * BOX_WIDTH) + ',' +
                        (genes.indexOf(d.gene) * BOX_HEIGHT + MAIN_PLOT_OFFSET) + ')';
                });

                // UPDATE rect (main)
                box.select('.sr-vm-expr-box')
                    .style('stroke', function(d) {
                        var zScore = parseFloat(d.zscore);
                        if (! zScore) {
                            return null;
                        } else if (Math.abs(zScore) > parseFloat(exprNum.value)) {
                            return zScore < 0 ? 'red' : 'green';
                        }
                        return null;
                    });
            }
            drawBoxes();

            (function drawLabels() {
                // DATA JOIN
                var geneLabel = svg.selectAll('.sr-vm-gene-label')
                    .data(genes);

                geneLabel.enter()
                    .append('text')
                    .attr('y', function(d) { return genes.indexOf(d) * BOX_HEIGHT + BOX_HEIGHT / 2 + MAIN_PLOT_OFFSET; })
                    .attr('x', - margin.left)
                    .attr('text-anchor', 'start')
                    .attr('font-size', BOX_HEIGHT / 2)
                    .text(function(d) { return d; });

                var subjectLabel = svg.selectAll('.sr-vm-subject-label')
                    .data(subjects);

                subjectLabel.enter()
                    .append('text')
                    .attr('transform', function(d) {
                        return 'translate(' + (subjects.indexOf(d) * BOX_WIDTH + BOX_WIDTH * 0.75) + ',' +
                            (height + 5) + ')rotate(-90)';
                    })
                    .attr('text-anchor', 'end')
                    .attr('font-size', BOX_WIDTH)
                    .text(function(d) { return d; });
            })();

            (function drawClinicalBoxes() {
                var catFeatures = clinicalFeatures.filter(function(d) { return d.indexOf('categoric') !== -1; });
                var numFeatures = clinicalFeatures.filter(function(d) { return d.indexOf('numeric') !== -1; });
                var catData = [];
                var numData = [];
                subjects.forEach(function(subject) {
                    bySubject.filterExact(subject);
                    var hit = bySubject.top(Infinity)[0];
                    clinicalFeatures.forEach(function(feature) {
                        if (hit && hit[feature]) {
                            if (catFeatures.indexOf(feature) !== -1) {
                                catData.push({
                                    subject: subject,
                                    feature: feature,
                                    value: hit[feature],
                                    type: 'categoric'
                                });
                            } else {
                                if (numFeatures.indexOf(feature) !== -1) {
                                    numData.push({
                                        subject: subject,
                                        feature: feature,
                                        value: hit[feature],
                                        type: 'numeric'
                                    });
                                }
                            }
                        }
                    });
                    bySubject.filterAll();
                });

                var scales = {};
                numFeatures.forEach(function(feature) {
                    var values = numData.filter(function(d) { return d.feature === feature; })
                        .map(function(d) { return parseInt(d.value); });
                    scales[feature] = d3.scale.linear()
                        .domain(d3.extent(values))
                        .range([0, BOX_WIDTH]);
                });

                // DATA JOIN
                var clinicalBox = svg.selectAll('.sr-vm-clinical-box')
                    .data(catData.concat(numData));

                // ENTER g
                var clinicalBoxEnter = clinicalBox.enter()
                    .append('g')
                    .attr('class', 'sr-vm-clinical-box');

                // ENTER rect
                clinicalBoxEnter.append('rect')
                    .attr('class', function(d) {
                        if (d.type === 'categoric') {
                            return 'sr-vm-cat-box';
                        } else {
                            return 'sr-vm-num-box';
                        }
                    })
                    .attr('height', CLINICAL_BOX_HEIGHT)
                    .attr('width', function(d) {
                        if (d.type === 'categoric') {
                            return BOX_WIDTH;
                        } else {
                            return scales[d.feature](d.value);
                        }
                    });

                // UPDATE g
                clinicalBox.attr('transform', function(d) {
                    return 'translate(' + (subjects.indexOf(d.subject) * BOX_WIDTH) + ',' +
                        (clinicalFeatures.indexOf(d.feature) * CLINICAL_BOX_HEIGHT) + ')';
                });
            });

            function getValuesForDimension(dimension, ascendingOrder) {
                var values = [];
                if (typeof ascendingOrder === 'undefined' || !ascendingOrder) {
                    values = dimension.top(Infinity).map(function(record) { return dimension.accessor(record); });
                } else {
                    values =dimension.bottom(Infinity).map(function(record) { return dimension.accessor(record); });
                }

                return values;
            }
        }

    }]);

