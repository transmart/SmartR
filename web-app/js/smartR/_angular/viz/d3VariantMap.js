//# sourceURL=d3VariantMap.js

'use strict';

window.smartRApp.directive('variantMap', [
    'smartRUtils',
    function(smartRUtils) {

        return {
            restrict: 'E',
            scope: {
                data: '=',
                width: '@',
                height: '@'
            },
            link: function (scope, element) {

                scope.$watch('data', function() {
                    $(element[0]).empty();
                    if (! $.isEmptyObject(scope.data)) {
                        smartRUtils.prepareWindowSize(scope.width, scope.height);
                        createVariantMap(scope, element[0]);
                    }
                });
            }
        };

        function createVariantMap(scope, viz) {
            var cf = crossfilter(scope.data.data);
            var bySubject = cf.dimension(function(d) { return d.subject; });
            var byGene = cf.dimension(function(d) { return d.gene; });

            var subjects = smartRUtils.unique(getValuesForDimension(bySubject));
            var genes = smartRUtils.unique(getValuesForDimension(byGene));

            var BOX_SIZE = 40;
            var height = genes.length * BOX_SIZE;
            var width = subjects.length * BOX_SIZE;

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

            (function drawGrid() {
                var subjects = smartRUtils.unique(getValuesForDimension(bySubject));
                var genes = smartRUtils.unique(getValuesForDimension(byGene));

                // DATA JOIN
                var verticalGridLine = svg.selectAll('.sr-vm-vgrid')
                    .data(subjects);

                // ENTER g
                var verticalGridLineEnter = verticalGridLine.enter()
                    .append('g')
                    .attr('class', 'sr-vm-vgrid sr-vm-grid');

                // ENTER line
                verticalGridLineEnter.append('line')
                    .attr('y1', 0)
                    .attr('y2', height);

                // ENTER text
                verticalGridLineEnter.append('text')
                    .attr('transform', 'translate(0,' + (height + 10) + ')rotate(45)')
                    .attr('text-anchor', 'start')
                    .style('font-size', BOX_SIZE / 2 + 'px')
                    .text(function(d) { return d; });

                // UPDATE g
                verticalGridLine.attr('transform', function(d) {
                    return 'translate(' + ((subjects.indexOf(d) + 0.5) * BOX_SIZE) + ', 0)';
                });

                // DATA JOIN
                var horizontalGridLine = svg.selectAll('.sr-vm-hgrid')
                    .data(genes);

                // ENTER g
                var horizontalGridLineEnter = horizontalGridLine.enter()
                    .append('g')
                    .attr('class', 'sr-vm-hgrid sr-vm-grid');

                // ENTER line
                horizontalGridLineEnter.append('line')
                    .attr('x1', 0)
                    .attr('x2', width);

                // ENTER text
                horizontalGridLineEnter.append('text')
                    .attr('transform', 'translate(-5, ' + (BOX_SIZE / 4) + ')')
                    .attr('text-anchor', 'end')
                    .style('font-size', (BOX_SIZE / 2) + 'px')
                    .text(function(d) { return d; });

                // UPDATE g
                horizontalGridLine.attr('transform', function(d) {
                    return 'translate(0, ' + ((genes.indexOf(d) + 0.5) * BOX_SIZE) + ')';
                });
            })();

            var redGreenScale = d3.scale.quantile()
                .domain([0, 1])
                .range(function() {
                    var colorSet = [];
                    var NUM = 100;
                    var i = NUM;
                    while (i--) {
                        colorSet.push(d3.rgb((255 * i) / NUM, 0, 0));
                    }
                    i = NUM;
                    while (i--) {
                        colorSet.push(d3.rgb(0, (255 * (NUM - i)) / NUM, 0));
                    }
                    return colorSet.reverse();
                }());

            (function drawBoxes() {
                var subjects = smartRUtils.unique(getValuesForDimension(bySubject));
                var genes = smartRUtils.unique(getValuesForDimension(byGene));
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
                        var aachange = !!consequences.filter(function(consequence) { return consequence !== 'synonymous'; }).length;
                        boxData.push({
                            subject: subject,
                            gene: gene,
                            zscore: zscore,
                            consequences: consequences,
                            aachange: aachange,
                            variants: hits.length,
                        });
                        byGene.filterAll();
                    });
                    bySubject.filterAll();
                });
                // DATA JOIN
                var box = svg.selectAll('.sr-vm-box')
                    .data(boxData, function(d) { return d.subject + ' ' + d.gene; });

                // ENTER g
                var boxEnter = box.enter()
                    .append('g')
                    .attr('class', 'sr-vm-box');

                // ENTER rect (main)
                boxEnter.append('rect')
                    .attr('class', 'sr-vm-main-box')
                    .attr('height', BOX_SIZE)
                    .attr('width', BOX_SIZE);

                // ENTER rect (expression box)
                boxEnter.append('rect')
                    .attr('class', 'sr-vm-expr-box')
                    .attr('height', BOX_SIZE / 2)
                    .attr('width', BOX_SIZE / 2)
                    .attr('x', BOX_SIZE / 2)
                    .style('fill', function(d) {
                        return redGreenScale(1 / (1 + Math.pow(Math.E, -d.zscore)));
                    });

                // ENTER line (consequence)
                boxEnter.append('line')
                    .attr('class', 'sr-vm-consq-line')
                    .attr('x1', function(d) { return d.aachange ? BOX_SIZE / 4 : 0; })
                    .attr('x2', function(d) { return d.aachange ? BOX_SIZE / 4 : BOX_SIZE / 2; })
                    .attr('y1', function(d) { return d.aachange ? BOX_SIZE / 2 : BOX_SIZE * 3 / 4; })
                    .attr('y2', function(d) { return d.aachange ? BOX_SIZE : BOX_SIZE * 3 / 4; })
                    .style('stroke', function(d) { return d.aachange ? '#F00' : '#0F0'; })
                    .style('visibility', function(d) { return d.variants ? 'visible' : 'hidden'; });

                // ENTER text (variants)
                boxEnter.append('text')
                    .attr('x', BOX_SIZE * 3 / 4)
                    .attr('y', BOX_SIZE - 3)
                    .style('font-size', BOX_SIZE / 2 + 'px')
                    .style('text-anchor', 'middle')
                    .text(function(d) { return d.variants ? d.variants : ''; });

                // UPDATE g
                box.attr('transform', function(d) {
                    return 'translate(' + (subjects.indexOf(d.subject) * BOX_SIZE) + ',' + (genes.indexOf(d.gene) * BOX_SIZE) + ')';
                });
            })();

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

