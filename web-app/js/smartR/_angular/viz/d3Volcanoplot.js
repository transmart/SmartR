//# sourceURL=d3Volcanoplot.js

'use strict';

window.smartRApp.directive('volcanoPlot', [
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
            templateUrl: $rootScope.smartRPath +  '/js/smartR/_angular/templates/volcanoplot.html',
            link: function (scope, element) {
                var viz = element.children()[1];
                scope.$watch('data', function() {
                    $(viz).empty();
                    if (!$.isEmptyObject(scope.data)) {
                        smartRUtils.prepareWindowSize(scope.width, scope.height);
                        createVolcanoplot(scope, viz);
                    }
                });
            }
        };

        function createVolcanoplot(scope, root) {
            var keggBtn = smartRUtils.getElementWithoutEventListeners('sr-volcanoplot-kegg-btn');
            keggBtn.addEventListener('click', function() {
                var labels = d3.selectAll('.point').filter(function(d) {
                    return d.negativeLog10PValues >= currentNegLog10P && Math.abs(d.logFC) >= currentLogFC;
                }).data().map(function(d) { return d.uid; });
                var genes = [];
                labels.forEach(function(label) {
                    var split = label.split("--");
                    split.shift();
                    genes = genes.concat(split);
                });

                var request = $.ajax({
                    url: pageInfo.basePath + '/SmartR/biocompendium',
                    type: 'POST',
                    timeout: 5000,
                    data: {
                        genes: genes.join(' ')
                    }
                });

                request.then(
                    function(response) {
                        var sessionID = response.match(/tmp_\d+/)[0];
                        var url = 'http://biocompendium.embl.de/' +
                            'cgi-bin/biocompendium.cgi?section=pathway&pos=0&background=whole_genome&session=' +
                            sessionID + '&list=gene_list_1__1&list_size=15&org=human';
                        window.open(url);
                    },
                    function(response) { alert("Error:", response); }
                );
            });

            var uids = scope.data.uids;
            var pValues = scope.data.pvalValues;
            var negativeLog10PValues = scope.data.negativeLog10PvalValues;
            var logFCs = scope.data.logfoldValues;

            var points = negativeLog10PValues.map(function (d, i) {
                return {
                    uid: uids[i],
                    pValue: pValues[i],
                    negativeLog10PValues: negativeLog10PValues[i],
                    logFC: logFCs[i]
                };
            });

            var currentLogFC = 0.5;
            var currentNegLog10P = -Math.log10(0.05);

            var margin = {top: 100, right: 100, bottom: 100, left: 100};
            var width = scope.width - margin.left - margin.right;
            var height = scope.height - margin.top - margin.bottom;

            var volcanoplot = d3.select(root).append('svg')
                .style('float', 'left')
                .attr('width', width + margin.left + margin.right)
                .attr('height', height + margin.top + margin.bottom)
                .append('g')
                .attr('transform', 'translate(' + margin.left + ',' + margin.top + ')');

            var tip = d3.tip()
                .attr('class', 'd3-tip')
                .offset([-10, 0])
                .html(function(d) { return d; });

            volcanoplot.call(tip);

            var maxAbsLogFCs = Math.max.apply(null, logFCs.map(Math.abs))

            var x = d3.scale.linear()
                .domain([-maxAbsLogFCs, maxAbsLogFCs])
                .range([0, width]);

            var y = d3.scale.linear()
                .domain(d3.extent(negativeLog10PValues))
                .range([height, 0]);

            var xAxis = d3.svg.axis()
                .scale(x)
                .orient('bottom');

            var yAxis = d3.svg.axis()
                .scale(y)
                .orient('left');

            volcanoplot.append('g')
                .attr('class', 'axis')
                .attr('transform', 'translate(0,' + height + ')')
                .call(xAxis);

            volcanoplot.append('g')
                .attr('class', 'axis')
                .call(yAxis);

            volcanoplot.append('g')
                .attr('class', 'x axis')
                .attr('transform', 'translate(0, 0)')
                .call(d3.svg.axis()
                    .scale(x)
                    .ticks(10)
                    .tickFormat('')
                    .innerTickSize(height)
                    .orient('bottom'));

            volcanoplot.append('g')
                .attr('class', 'y axis')
                .attr('transform', 'translate(' + width + ',' + 0 + ')')
                .call(d3.svg.axis()
                    .scale(y)
                    .ticks(10)
                    .tickFormat('')
                    .innerTickSize(width)
                    .orient('left'));

            volcanoplot.append('text')
                .attr('class', 'text axisText')
                .attr('x', width / 2)
                .attr('y', height + 40)
                .attr('text-anchor', 'middle')
                .text('log2 FC');

            volcanoplot.append('text')
                .attr('class', 'text axisText')
                .attr('text-anchor', 'middle')
                .attr('transform', 'translate(' + (-40) + ',' + (height / 2) + ')rotate(-90)')
                .text('- log10 p');

            function pDragged() {
                var yPos = d3.event.y;
                if (yPos < 0) {
                    yPos = 0;
                }
                if (yPos > height) {
                    yPos = height;
                }

                d3.selectAll('.pLine')
                    .attr('y1', yPos)
                    .attr('y2', yPos);

                d3.selectAll('.pHandle')
                    .attr('y', yPos - 6);

                d3.selectAll('.pText')
                    .attr('y', yPos)
                    .text('p = ' + (1 / Math.pow(10, y.invert(yPos))).toFixed(5));

                currentNegLog10P = y.invert(yPos);

                d3.selectAll('.point')
                    .style('fill', function (d) {
                        return getColor(d);
                    });

                drawVolcanotable(getTopRankedPoints().data());
            }

            var pDrag = d3.behavior.drag()
                .on('drag', pDragged)
                .on('dragend', d3.tip);

            volcanoplot.append('line')
                .attr('class', 'pLine')
                .attr('x1', 0)
                .attr('y1', y(currentNegLog10P))
                .attr('x2', width)
                .attr('y2', y(currentNegLog10P));

            volcanoplot.append('rect')
                .attr('class', 'pHandle')
                .attr('x', 0)
                .attr('y', y(currentNegLog10P) - 6)
                .attr('width', width)
                .attr('height', 12)
                .call(pDrag);

            volcanoplot.append('text')
                .attr('class', 'text pText')
                .attr('x', width + 5)
                .attr('y', y(currentNegLog10P))
                .attr('dy', '0.35em')
                .attr('text-anchor', 'start')
                .text('p = 0.0500')
                .style('fill', 'red');

            function lFCDragged() {
                var xPos = d3.event.x;

                if (xPos < 0) {
                    xPos = 0;
                }
                if (xPos > width) {
                    xPos = width;
                }

                var logFC = x.invert(xPos);

                d3.selectAll('.logFCLine.left')
                    .attr('x1', x(-Math.abs(logFC)))
                    .attr('x2', x(-Math.abs(logFC)));
                d3.selectAll('.logFCHandle.left')
                    .attr('x', x(-Math.abs(logFC)) - 6);
                d3.selectAll('.logFCText.left')
                    .attr('x', x(-Math.abs(logFC)))
                    .text('log2FC = ' + (-Math.abs(logFC)).toFixed(2));

                d3.selectAll('.logFCLine.right')
                    .attr('x1', x(Math.abs(logFC)))
                    .attr('x2', x(Math.abs(logFC)));
                d3.selectAll('.logFCHandle.right')
                    .attr('x', x(Math.abs(logFC)) - 6);
                d3.selectAll('.logFCText.right')
                    .attr('x', x(Math.abs(logFC)))
                    .text('log2FC = ' + Math.abs(logFC).toFixed(2));

                currentLogFC = Math.abs(logFC);

                d3.selectAll('.point')
                    .style('fill', function (d) {
                        return getColor(d);
                    });

                drawVolcanotable(getTopRankedPoints().data());
            }

            var lFCDrag = d3.behavior.drag()
                .on('drag', lFCDragged);

            volcanoplot.append('line')
                .attr('class', 'left logFCLine')
                .attr('x1', x(-currentLogFC))
                .attr('y1', height)
                .attr('x2', x(-currentLogFC))
                .attr('y2', 0);

            volcanoplot.append('rect')
                .attr('class', 'left logFCHandle')
                .attr('x', x(-currentLogFC) - 6)
                .attr('y', 0)
                .attr('width', 12)
                .attr('height', height)
                .call(lFCDrag);

            volcanoplot.append('text')
                .attr('class', 'text left logFCText')
                .attr('x', x(-currentLogFC))
                .attr('y', -15)
                .attr('dy', '0.35em')
                .attr('text-anchor', 'end')
                .text('log2FC = ' + (-currentLogFC))
                .style('fill', '#0000FF');

            volcanoplot.append('line')
                .attr('class', 'right logFCLine')
                .attr('x1', x(currentLogFC))
                .attr('y1', height)
                .attr('x2', x(currentLogFC))
                .attr('y2', 0);

            volcanoplot.append('rect')
                .attr('class', 'right logFCHandle')
                .attr('x', x(currentLogFC) - 6)
                .attr('y', 0)
                .attr('width', 12)
                .attr('height', height)
                .call(lFCDrag);

            volcanoplot.append('text')
                .attr('class', 'text right logFCText')
                .attr('x', x(currentLogFC))
                .attr('y', -15)
                .attr('dy', '0.35em')
                .attr('text-anchor', 'start')
                .text('log2FC = ' + currentLogFC)
                .style('fill', '#0000FF');

            function getTopRankedPoints() {
                return d3.selectAll('.point').filter(function (d) {
                    return d.negativeLog10PValues > currentNegLog10P && Math.abs(d.logFC) > currentLogFC;
                });
            }

            function getColor(point) {
                if (point.negativeLog10PValues < currentNegLog10P && Math.abs(point.logFC) < currentLogFC) {
                    return '#000000';
                }
                if (point.negativeLog10PValues >= currentNegLog10P && Math.abs(point.logFC) < currentLogFC) {
                    return '#FF0000';
                }
                if (point.negativeLog10PValues >= currentNegLog10P && Math.abs(point.logFC) >= currentLogFC) {
                    return '#00FF00';
                }
                return '#0000FF';
            }

            function resetVolcanotable() {
                d3.selectAll('.volcanoplot-table').remove();
            }

            function drawVolcanotable(points) {
                resetVolcanotable();
                if (!points.length) {
                    return;
                }
                var columns = ['uid', 'logFC', 'negativeLog10PValues', 'pValue'];
                var HEADER = ['ID', 'log2 FC', '- log10 p', 'p'];
                var table = d3.select(root).append('table')
                    .attr('class', 'volcanoplot-table');
                var thead = table.append('thead');
                var tbody = table.append('tbody');

                thead.append('tr')
                    .attr('class', 'mytr')
                    .selectAll('th')
                    .data(HEADER)
                    .enter()
                    .append('th')
                    .attr('class', 'myth')
                    .text(function (d) {
                        return d;
                    });

                var rows = tbody.selectAll('tr')
                    .data(points)
                    .enter()
                    .append('tr')
                    .attr('class', 'mytr');

                rows.selectAll('td')
                    .data(function (row) {
                        return columns.map(function (column) {
                            return {column: column, value: row[column]};
                        });
                    })
                    .enter()
                    .append('td')
                    .attr('class', 'text mytd')
                    .text(function (d) {
                        return d.value;
                    });
            }

            function updateVolcano() {
                var point = volcanoplot.selectAll('.point')
                    .data(points, function (d) {
                        return d.uid;
                    });

                point.enter()
                    .append('rect')
                    .attr('class', function(d) { return 'point uid-' + smartRUtils.makeSafeForCSS(d.uid); })
                    .attr('x', function (d) {
                        return x(d.logFC) - 2;
                    })
                    .attr('y', function (d) {
                        return y(d.negativeLog10PValues) - 2;
                    })
                    .attr('width', 4)
                    .attr('height', 4)
                    .style('fill', function (d) {
                        return getColor(d);
                    })
                    .on('mouseover', function (d) {
                        var html = 'ID: ' + d.uid + '<br/>' +
                            'p-value: ' + d.pValue + '<br/>' +
                            '-log10 p: ' + d.negativeLog10PValues + '<br/>' +
                            'log2FC: ' + d.logFC;
                        tip.show(html);
                    })
                    .on('mouseout', function () {
                        tip.hide();
                    });

                point.exit()
                    .attr('r', 0)
                    .remove();
            }

            updateVolcano();
            drawVolcanotable(getTopRankedPoints().data());
        }

    }]);

