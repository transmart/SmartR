//# sourceURL=d3Correlation.js

'use strict';

window.smartRApp.directive('correlationPlot', [
    'smartRUtils',
    'rServeService',
    function(smartRUtils, rServeService) {

        return {
            restrict: 'E',
            scope: {
                data: '=',
                width: '@',
                height: '@'
            },
            link: function (scope, element) {

                /**
                 * Watch data model (which is only changed by ajax calls when we want to (re)draw everything)
                 */
                scope.$watch('data', function() {
                    $(element[0]).empty();
                    if (! $.isEmptyObject(scope.data)) {
                        smartRUtils.prepareWindowSize(scope.width, scope.height);
                        createCorrelationViz(scope, element[0]);
                    }
                });
            }
        };

        function createCorrelationViz(scope, root) {
            var animationDuration = 500;
            var bins = 10;
            var w = parseInt(scope.width);
            var h = parseInt(scope.height);
            var margin = {top: 20, right: 250, bottom: h / 4 + 230, left: w / 4};
            var width = w * 3 / 4 - margin.left - margin.right;
            var height = h * 3 / 4 - margin.top - margin.bottom;
            var bottomHistHeight = margin.bottom;
            var leftHistHeight = margin.left;
            var colors = ['#33FF33', '#3399FF', '#CC9900', '#CC99FF', '#FFFF00', 'blue'];
            var x = d3.scale.linear()
                .domain(d3.extent(scope.data.points, function(d) { return d.x; }))
                .range([0, width]);
            var y = d3.scale.linear()
                .domain(d3.extent(scope.data.points, function(d) { return d.y; }))
                .range([height, 0]);

            var annotations = scope.data.annotations.sort();
            for (var i = 0; i < annotations.length; i++) {
                var annotation = annotations[i];
                if (annotation === '') {
                    annotations[i] = 'Default';
                }
            }

            var xArrLabel = scope.data.xArrLabel[0];
            var yArrLabel = scope.data.yArrLabel[0];

            var correlation,
                pvalue,
                regLineSlope,
                regLineYIntercept,
                patientIDs,
                points,
                method,
                transformation,
                minX,
                maxX,
                minY,
                maxY;
            function setData(data) {
                correlation = data.correlation[0];
                pvalue = data.pvalue[0];
                regLineSlope = data.regLineSlope[0];
                regLineYIntercept = data.regLineYIntercept[0];
                method = data.method[0];
                transformation = data.transformation[0];
                patientIDs = data.patientIDs;
                points = data.points;
                var xValues =  data.points.map(function(d) { return d.x; });
                var yValues =  data.points.map(function(d) { return d.y; });
                minX = Math.min.apply(null, xValues);
                minY = Math.min.apply(null, yValues);
                maxX = Math.max.apply(null, xValues);
                maxY = Math.max.apply(null, yValues);
            }

            setData(scope.data);

            function updateStatistics(patientIDs, scatterUpdate, init) {
                if (! init) {
                    patientIDs = patientIDs.length !== 0 ? patientIDs : d3.selectAll('.point').map(function(d) {
                        return d.patientID;
                    });
                }
                var args = { method: method, transformation: transformation, selectedPatientIDs: patientIDs };

                rServeService.startScriptExecution({
                    taskType: 'run',
                    arguments: args
                }).then(
                    function (response) {
                        var results = JSON.parse(response.result.artifacts.value);
                        if (init) {
                            removePlot();
                            scope.data = results;
                        } else {
                            setData(results);
                            if (scatterUpdate) {
                                updateScatterplot();
                            }
                            updateRegressionLine();
                            updateLegends();
                            updateHistogram();
                        }
                    },
                    function (response) {
                        console.error(response);
                    }
                );
            }

            function removePlot() {
                d3.select(root).selectAll('*').remove();
                d3.selectAll('.d3-tip').remove();
            }

            var svg = d3.select(root).append('svg')
                .attr('width', width + margin.left + margin.right)
                .attr('height', height + margin.top + margin.bottom)
                .append('g')
                .attr('transform', 'translate(' + margin.left + ',' + margin.top + ')')
                .on('contextmenu', function() {
                    d3.event.preventDefault();
                    contextMenu.show('<input id="sr-correlation-zoom-btn" value="Zoom" class="sr-ctx-menu-btn"><br/>' +
                                     '<input id="sr-correlation-exclude-btn" value="Exclude" class="sr-ctx-menu-btn"><br/>' +
                                     '<input id="sr-correlation-reset-btn" value="Reset" class="sr-ctx-menu-btn">');
                    _addEventListeners();
                });

            function _addEventListeners() {
                smartRUtils.getElementWithoutEventListeners('sr-correlation-zoom-btn').addEventListener('click', function() {
                    contextMenu.hide();
                    zoomSelection();
                });
                smartRUtils.getElementWithoutEventListeners('sr-correlation-exclude-btn').addEventListener('click', function() {
                    contextMenu.hide();
                    excludeSelection();
                });
                smartRUtils.getElementWithoutEventListeners('sr-correlation-reset-btn').addEventListener('click', function() {
                    contextMenu.hide();
                    reset();
                });
            }

            var tip = d3.tip()
                .attr('class', 'd3-tip')
                .offset([-10, 0])
                .html(function(d) { return d; });

            svg.call(tip);

            svg.append('g')
                .attr('class', 'x axis')
                .attr('transform', 'translate(0, 0)')
                .call(d3.svg.axis()
                    .scale(x)
                    .ticks(10)
                    .tickFormat('')
                    .innerTickSize(height)
                    .orient('bottom'));

            svg.append('text')
                .attr('class', 'axisLabels')
                .attr('transform', 'translate(' + width / 2 + ',' + (- 10) + ')')
                .text(smartRUtils.shortenConcept(xArrLabel) + ' (' + transformation + ')');

            svg.append('g')
                .attr('class', 'y axis')
                .attr('transform', 'translate(' + width + ',' + 0 + ')')
                .call(d3.svg.axis()
                    .scale(y)
                    .ticks(10)
                    .tickFormat('')
                    .innerTickSize(width)
                    .orient('left'));

            svg.append('text')
                .attr('class', 'axisLabels')
                .attr('transform', 'translate('  + (width + 10) + ',' + height / 2 + ')rotate(90)')
                .text(smartRUtils.shortenConcept(yArrLabel) + ' (' + transformation + ')');

            svg.append('g')
                .attr('class', 'x axis')
                .attr('transform', 'translate(' + 0 + ',' + height + ')')
                .call(d3.svg.axis()
                    .scale(x)
                    .orient('top'));

            svg.append('g')
                .attr('class', 'y axis')
                .attr('transform', 'translate(' + 0 + ',' + 0 + ')')
                .call(d3.svg.axis()
                    .scale(y)
                    .orient('right'));

            function excludeSelection() {
                var remainingPatientIDs = d3.selectAll('.point:not(.selected)').map(function(d) {
                    return d.patientID;
                });
                updateStatistics(remainingPatientIDs, true);
            }

            function zoomSelection() {
                if (d3.selectAll('.point.selected').size() < 2) {
                    alert('Please select at least two elements before zooming!');
                    return;
                }
                var selectedPatientIDs = d3.selectAll('.point.selected').map(function(d) { return d.patientID; });
                updateStatistics(selectedPatientIDs, false, true);
            }

            var contextMenu = d3.tip()
                .attr('class', 'd3-tip sr-contextmenu')
                .offset([-10, 0])
                .html(function(d) { return d; });

            svg.call(contextMenu);

            function updateSelection() {
                var extent = brush.extent();
                var x0 = x.invert(extent[0][0]);
                var x1 = x.invert(extent[1][0]);
                var y0 = y.invert(extent[0][1]);
                var y1 = y.invert(extent[1][1]);
                svg.selectAll('.point')
                    .classed('selected', false)
                    .style('fill', function(d) { return getColor(d.annotation); })
                    .style('stroke', 'white')
                    .filter(function(d) {
                        return x0 <= d.x && d.x <= x1 && y1 <= d.y && d.y <= y0;
                    })
                    .classed('selected', true)
                    .style('fill', 'white')
                    .style('stroke', function(d) { return getColor(d.annotation); });
            }

            var brush = d3.svg.brush()
                .x(d3.scale.identity().domain([0, width]))
                .y(d3.scale.identity().domain([0, height]))
                .on('brushend', function() {
                    contextMenu.hide();
                    updateSelection();
                    var selectedPatientIDs = d3.selectAll('.point.selected').map(function(d) { return d.patientID; });
                    updateStatistics(selectedPatientIDs);
                });

            svg.append('g')
                .attr('class', 'brush')
                .on('mousedown', function() {
                    return d3.event.button === 2 ? d3.event.stopImmediatePropagation() : null;
                })
                .call(brush);

            function getColor(annotation) {
                return annotation && annotation !== 'Default' ? colors[annotations.indexOf(annotation)] : 'black';
            }

            function updateScatterplot() {
                var point = svg.selectAll('.point')
                    .data(points, function(d) { return d.patientID; });

                point.enter()
                    .append('circle')
                    .attr('class', 'point')
                    .attr('cx', function(d) { return x(d.x); })
                    .attr('cy', function(d) { return y(d.y); })
                    .attr('r', 5)
                    .style('fill', function(d) { return getColor(d.annotation); })
                    .on('mouseover', function(d) {
                        d3.select(this).style('fill', '#FF0000');
                        tip.show(smartRUtils.shortenConcept(xArrLabel) + ': ' + d.x + '<br/>' +
                            smartRUtils.shortenConcept(yArrLabel) + ': ' + d.y + '<br/>' +
                            'Patient ID: ' + d.patientID + '<br/>' +
                            (d.annotation ? 'Tag: ' + d.annotation : ''));
                    })
                    .on('mouseout', function() {
                        var p = d3.select(this);
                        p.style('fill', function(d) {
                            return p.classed('selected') ? '#FFFFFF' : getColor(d.annotation);
                        });
                        tip.hide();
                    });

                point.exit()
                    .classed('selected', false)
                    .transition()
                    .duration(animationDuration)
                    .attr('r', 0)
                    .remove();
            }

            function updateHistogram() {
                var bottomHistData = d3.layout.histogram()
                    .bins(bins)(points.map(function(d) { return d.x; }));
                var leftHistData = d3.layout.histogram()
                    .bins(bins)(points.map(function(d) { return d.y; }));

                var bottomHistHeightScale = d3.scale.linear()
                    .domain([0, Math.max.apply(null, bottomHistData.map(function(d) { return d.y; }))])
                    .range([1, bottomHistHeight]);
                var leftHistHeightScale = d3.scale.linear()
                    .domain([0, Math.max.apply(null, leftHistData.map(function(d) { return d.y; }))])
                    .range([2, leftHistHeight]);

                var bottomHistGroup = svg.selectAll('.bar.bottom')
                    .data(Array(bins).fill().map(function(_, i) { return i; }));
                var bottomHistGroupEnter = bottomHistGroup.enter()
                    .append('g')
                    .attr('class', 'bar bottom');
                var bottomHistGroupExit = bottomHistGroup.exit();

                bottomHistGroupEnter.append('rect')
                    .attr('y', height + 1);
                bottomHistGroup.selectAll('rect')
                    .transition()
                    .delay(function(d) { return d * 25; })
                    .duration(animationDuration)
                    .attr('x', function(d) { return x(bottomHistData[d].x); })
                    .attr('width', function() { return (x(maxX) - x(minX)) / bins; })
                    .attr('height', function(d) { return bottomHistHeightScale(bottomHistData[d].y) - 1; });
                bottomHistGroupExit.selectAll('rect')
                    .transition()
                    .duration(animationDuration)
                    .attr('height', 0);

                bottomHistGroupEnter.append('text')
                    .attr('dy', '.35em')
                    .attr('text-anchor', 'middle');
                bottomHistGroup.selectAll('text')
                    .text(function(d) { return bottomHistData[d].y || ''; })
                    .transition()
                    .delay(function(d) { return d * 25; })
                    .duration(animationDuration)
                    .attr('x', function(d) { return x(bottomHistData[d].x) + (x(maxX) - x(minX)) / bins / 2; })
                    .attr('y', function(d) { return height + bottomHistHeightScale(bottomHistData[d].y) - 10; });
                bottomHistGroupExit.selectAll('text')
                    .text('');

                var leftHistGroup = svg.selectAll('.bar.left')
                    .data(Array(bins).fill().map(function(_, i) { return i; }));
                var leftHistGroupEnter = leftHistGroup.enter()
                    .append('g')
                    .attr('class', 'bar left');
                var leftHistGroupExit = leftHistGroup.exit();

                leftHistGroupEnter.append('rect');
                leftHistGroup.selectAll('rect')
                    .transition()
                    .delay(function(d) { return d * 25; })
                    .duration(animationDuration)
                    .attr('x', function(d) { return - leftHistHeightScale(leftHistData[d].y) + 1; })
                    .attr('y', function(d) { return y(leftHistData[d].x) - (y(minY) - y(maxY))/ bins; })
                    .attr('width', function(d) { return leftHistHeightScale(leftHistData[d].y) - 2; })
                    .attr('height', function() { return (y(minY) - y(maxY))/ bins; });
                leftHistGroupExit.selectAll('rect')
                    .transition()
                    .duration(animationDuration)
                    .attr('height', 0);

                leftHistGroupEnter.append('text')
                    .attr('dy', '.35em')
                    .attr('text-anchor', 'middle');
                leftHistGroup.selectAll('text')
                    .text(function(d) { return leftHistData[d].y || ''; })
                    .transition()
                    .delay(function(d) { return d * 25; })
                    .duration(animationDuration)
                    .attr('x', function(d) { return - leftHistHeightScale(leftHistData[d].y) + 10; })
                    .attr('y', function(d) { return y(leftHistData[d].x) - (y(minY) - y(maxY))/ bins / 2; });
                leftHistGroupExit.selectAll('text')
                    .text('');
            }

            function updateLegends() {
                var _LEGEND_LINE_SPACING = 15,
                    _LEGEND_RECT_SIZE = 10;
                var text = [
                    {name: 'Correlation Coefficient: ' + correlation},
                    {name: 'p-value: ' + pvalue},
                    {name: 'Method: ' + method},
                    {name: 'Selected: ' + d3.selectAll('.point.selected').size() || d3.selectAll('.point').size()},
                    {name: 'Displayed: ' + d3.selectAll('.point').size()}
                ];

                annotations.forEach(function(annotation) {
                    text.push({color: getColor(annotation), name: annotation});
                });

                var legend = svg.selectAll('.legend')
                    .data(text, function(d) { return d.name; });

                var legendEnter = legend.enter()
                    .append('g')
                    .attr('class', 'legend')
                    .attr('transform', function(d, i) {
                        return 'translate(' + (width + 40) + ',' + (i * (_LEGEND_RECT_SIZE + _LEGEND_LINE_SPACING)) + ')';
                    });

                legendEnter.append('text');
                legendEnter.append('rect');

                legend.select('text')
                    .attr('x', function(d) { return d.color ? 20 : 0; })
                    .attr('y', 9)
                    .text(function(d) { return d.name; });

                legend.select('rect')
                    .attr('width', _LEGEND_RECT_SIZE)
                    .attr('height', _LEGEND_RECT_SIZE)
                    .style('fill', function(d) { return d.color; })
                    .style('visibility', function(d) { return d.color ? 'visible' : 'hidden'; });

                legend.exit()
                    .remove();
            }

            function updateRegressionLine() {
                var regressionLine = svg.selectAll('.regressionLine')
                    .data(regLineSlope === 'NA' ? [] : [0], function(d) { return d; });
                regressionLine.enter()
                    .append('line')
                    .attr('class', 'regressionLine')
                    .on('mouseover', function () {
                        d3.select(this).attr('stroke', 'red');
                        tip.show('slope: ' + regLineSlope + '<br/>intercept: ' + regLineYIntercept);
                    })
                    .on('mouseout', function () {
                        d3.select(this).attr('stroke', 'orange');
                        tip.hide();
                    });

                var x1 = x(minX),
                    y1 = y(regLineYIntercept + regLineSlope * minX),
                    x2 = x(maxX),
                    y2 = y(regLineYIntercept + regLineSlope * maxX);

                x1 = x1 < 0 ? 0 : x1;
                x1 = x1 > width ? width : x1;

                x2 = x2 < 0 ? 0 : x2;
                x2 = x2 > width ? width : x2;

                y1 = y1 < 0 ? 0 : y1;
                y1 = y1 > height ? height : y1;

                y2 = y2 < 0 ? 0 : y2;
                y2 = y2 > height ? height : y2;

                regressionLine.transition()
                    .duration(animationDuration)
                    .attr('x1', x1)
                    .attr('y1', y1)
                    .attr('x2', x2)
                    .attr('y2', y2);

                regressionLine.exit()
                    .remove();
            }

            function reset() {
                updateStatistics([], true, true);
            }

            updateScatterplot();
            updateHistogram();
            updateRegressionLine();
            updateLegends();
        }

    }]);
