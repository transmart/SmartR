//# sourceURL=d3Boxplot.js

'use strict';

window.smartRApp.directive('boxplot', [
    'smartRUtils',
    'rServeService',
    '$rootScope',
    function(smartRUtils, rServeService, $rootScope) {

    return {
        restrict: 'E',
        scope: {
            data: '=',
            width: '@',
            height: '@'
        },
        templateUrl: $rootScope.smartRPath +  '/js/smartR/_angular/templates/boxplot.html',
        link: function (scope, element) {
            var template_ctrl = element.children()[0],
                template_viz = element.children()[1];
            /**
             * Watch data model (which is only changed by ajax calls when we want to (re)draw everything)
             */
            scope.$watch('data', function () {
                $(template_viz).empty();
                if (! $.isEmptyObject(scope.data)) {
                    smartRUtils.prepareWindowSize(scope.width, scope.height);
                    scope.showControls = true;
                    createBoxplot(scope, template_viz, template_ctrl);
                }
            });
        }
    };

    function createBoxplot(scope, root) {
        var concept = '',
            globalMin = Number.MIN_VALUE,
            globalMax = Number.MAX_VALUE,
            categories = [],
            excludedPatientIDs = [],
            useLog = false;
        function setData(data) {
            concept = data.concept[0];
            globalMin = data.globalMin[0];
            globalMax = data.globalMax[0];
            categories = data['Subset 2'] ? ['Subset 1', 'Subset 2'] : ['Subset 1'];
            excludedPatientIDs = data.excludedPatientIDs;
            useLog = data.useLog[0];
        }
        setData(scope.data);

        var removeBtn = smartRUtils.getElementWithoutEventListeners('sr-boxplot-remove-btn');
        removeBtn.addEventListener('click', removeOutliers);

        var resetBtn = smartRUtils.getElementWithoutEventListeners('sr-boxplot-reset-btn');
        resetBtn.addEventListener('click', reset);

        var kdeCheck = smartRUtils.getElementWithoutEventListeners('sr-boxplot-kde-check');
        kdeCheck.addEventListener('change', function() { swapKDE(kdeCheck.checked); });
        kdeCheck.checked = false;

        var jitterCheck = smartRUtils.getElementWithoutEventListeners('sr-boxplot-jitter-check');
        jitterCheck.addEventListener('change', function() { swapJitter(jitterCheck.checked); });
        jitterCheck.checked = false;

        var logCheck = smartRUtils.getElementWithoutEventListeners('sr-boxplot-log-check');
        logCheck.checked = useLog;
        logCheck.addEventListener('change', function() { swapLog(logCheck.checked); });

        var animationDuration = 1000;

        var width = parseInt(scope.width);
        var height = parseInt(scope.height);
        var margin = {top: 20, right: 60, bottom: 200, left: 280};

        var boxWidth = 0.12 * width;
        var whiskerLength = boxWidth / 6;

        var colorScale = d3.scale.quantile()
            .range(['rgb(158,1,66)', 'rgb(213,62,79)', 'rgb(244,109,67)', 'rgb(253,174,97)', 'rgb(254,224,139)', 'rgb(255,255,191)', 'rgb(230,245,152)', 'rgb(171,221,164)', 'rgb(102,194,165)', 'rgb(50,136,189)', 'rgb(94,79,162)']);

        var boxplot = d3.select(root).append('svg')
            .attr('width', width + margin.left + margin.right)
            .attr('height', height + margin.top + margin.bottom)
            .append('g')
            .attr('transform', 'translate(' + margin.left + ',' + margin.top + ')');

        var x = d3.scale.ordinal()
            .domain(categories)
            .rangeBands([0, width], 1, 0.5);

        var y = d3.scale.linear()
            .domain([globalMin, globalMax])
            .range([height, 0]);

        var yAxis = d3.svg.axis()
            .scale(y)
            .orient('left');

        var xAxis = d3.svg.axis()
            .scale(x)
            .orient('bottom');

        boxplot.append('g')
            .attr('class', 'y axis')
            .attr('transform', 'translate(' + 0 + ',' + 0 + ')')
            .call(yAxis);

        boxplot.append('g')
            .attr('class', 'x axis')
            .attr('transform', 'translate(' + 0 + ',' + (height + 20) + ')')
            .call(xAxis)
            .selectAll('text')
            .attr('dy', '.35em')
            .attr('transform', 'translate(' + 0 + ',' + 5 + ')rotate(45)')
            .style('text-anchor', 'start');

        boxplot.append('text')
            .attr('transform', 'translate(' + (-40) + ',' + (height / 2) + ')rotate(-90)')
            .attr('text-anchor', 'middle')
            .text(smartRUtils.shortenConcept(concept) + (useLog ? ' (log2)' : ''));

        var tip = d3.tip()
            .attr('class', 'd3-tip')
            .offset([-10, 0])
            .html(function(d) { return d; });

        boxplot.call(tip);

        var brush = d3.svg.brush()
            .x(d3.scale.identity().domain([0, width]))
            .y(d3.scale.identity().domain([-5, height + 5]))
            .on('brush', function () {
                contextMenu.hide();
                updateSelection();
            });

        boxplot.append('g')
            .attr('class', 'brush')
            .on('mousedown', function () {
                if (d3.event.button === 2) {
                    d3.event.stopImmediatePropagation();
                }
            })
            .call(brush);

        var contextMenu = d3.tip()
            .attr('class', 'd3-tip sr-contextmenu')
            .html(function(d) { return d; });

        boxplot.call(contextMenu);

        function _addEventListeners() {
            smartRUtils.getElementWithoutEventListeners('sr-boxplot-exclude-btn').addEventListener('click', function() {
                contextMenu.hide();
                excludeSelection();
            });

            smartRUtils.getElementWithoutEventListeners('sr-boxplot-reset-btn-2').addEventListener('click', function() {
                contextMenu.hide();
                reset();
            });
        }

        boxplot.on('contextmenu', function () {
            d3.event.preventDefault();
            contextMenu.show('<input id="sr-boxplot-exclude-btn" value="Exclude" class="sr-ctx-menu-btn"><br/>' +
                             '<input id="sr-boxplot-reset-btn-2" value="Reset" class="sr-ctx-menu-btn">');
            _addEventListeners();
        });

        var currentSelection;
        function updateSelection() {
            d3.selectAll('.point').classed('brushed', false);
            var extent = brush.extent();
            var left = extent[0][0],
                top = extent[0][1],
                right = extent[1][0],
                bottom = extent[1][1];
            currentSelection = d3.selectAll('.point')
                .filter(function (d) {
                    var point = d3.select(this);
                    return y(d.value) >= top && y(d.value) <= bottom && point.attr('cx') >= left && point.attr('cx') <= right;
                })
                .classed('brushed', true)
                .map(function(d) { return d.patientID; });
        }

        function excludeSelection() {
            excludedPatientIDs = excludedPatientIDs.concat(currentSelection);
            var settings = { excludedPatientIDs: excludedPatientIDs, useLog: logCheck.checked };

            rServeService.startScriptExecution({
                taskType: 'run',
                arguments: settings
            }).then(
                function (response) {
                    removePlot();
                    scope.data = JSON.parse(response.result.artifacts.value);
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

        function removeOutliers() {
            currentSelection = d3.selectAll('.outlier').map(function (d) { return d.patientID; });
            if (currentSelection) { excludeSelection(); }
        }

        function kernelDensityEstimator(kernel, x) {
            return function (sample) {
                return x.map(function (x) {
                    return [x, d3.mean(sample, function (v) {
                        return kernel(x - v);
                    })];
                });
            };
        }

        function epanechnikovKernel(scale) {
            return function (u) {
                return Math.abs(u /= scale) <= 1 ? 0.75 * (1 - u * u) / scale : 0;
            };
        }

        // function gaussKernel(scale) {
        //     return function (u) {
        //         return Math.exp(-u * u / 2) / Math.sqrt(2 * Math.PI) / scale;
        //     };
        // }

        function swapKDE(checked) {
            if (!checked) {
                d3.selectAll('.line').attr('visibility', 'hidden');
            } else {
                d3.selectAll('.line').attr('visibility', 'visible');
            }
        }

        function shortenNodeLabel(label) {
            label = label.replace(/\W+/g, '');
            return label;
        }

        var jitterWidth = 1.0;

        function swapJitter() {
            categories.forEach(function(category) {
                d3.selectAll('.point.' + smartRUtils.makeSafeForCSS(shortenNodeLabel(category)))
                    .transition()
                    .duration(animationDuration)
                    .attr('cx', function(d) {
                        return jitterCheck.checked ? x(category) + boxWidth * jitterWidth * d.jitter : x(category);
                    });
            });
        }

        function swapLog() {
            var settings = { excludedPatientIDs: excludedPatientIDs, useLog: logCheck.checked };
            rServeService.startScriptExecution({
                taskType: 'run',
                arguments: settings
            }).then(
                function (response) {
                    removePlot();
                    scope.data = JSON.parse(response.result.artifacts.value);
                },
                function (response) {
                    console.error(response);
                }
            );
        }

        var boxes = {};
        categories.forEach(function(category) {
            boxes[category] = boxplot.append('g');
            var params = scope.data[category];
            createBox(params, category, boxes[category]);
        });

        function createBox(params, category, box) {
            var boxLabel = shortenNodeLabel(category);
            colorScale.domain(d3.extent(y.domain()));

            var whisker = box.selectAll('.whisker')
                .data([params.upperWhisker, params.lowerWhisker], function (d, i) {
                    return boxLabel + '-whisker-' + i;
                });

            whisker.enter()
                .append('line')
                .attr('class', 'whisker');

            whisker.transition()
                .duration(animationDuration)
                .attr('x1', x(category) - whiskerLength / 2)
                .attr('y1', function (d) { return y(d); })
                .attr('x2', x(category) + whiskerLength / 2)
                .attr('y2', function (d) { return y(d); });

            var whiskerLabel = box.selectAll('.whiskerLabel')
                .data([params.upperWhisker, params.lowerWhisker], function (d, i) {
                    return boxLabel + '-whiskerLabel-' + i;
                });

            whiskerLabel.enter()
                .append('text')
                .attr('class', 'whiskerLabel boxplotValue');

            whiskerLabel.transition()
                .duration(animationDuration)
                .attr('x', x(category) + whiskerLength / 2)
                .attr('y', function (d) { return y(d); })
                .attr('dx', '.35em')
                .attr('dy', '.35em')
                .attr('text-anchor', 'start')
                .text(function (d) { return d; });

            var hingeLength = boxWidth;
            var hinge = box.selectAll('.hinge')
                .data([params.upperHinge, params.lowerHinge], function (d, i) {
                    return boxLabel + '-hinge-' + i;
                });

            hinge.enter()
                .append('line')
                .attr('class', 'hinge');

            hinge.transition()
                .duration(animationDuration)
                .attr('x1', x(category) - hingeLength / 2)
                .attr('y1', function (d) { return y(d); })
                .attr('x2', x(category) + hingeLength / 2)
                .attr('y2', function (d) { return y(d); });

            var hingeLabel = box.selectAll('.hingeLabel')
                .data([params.upperHinge, params.lowerHinge], function (d, i) {
                    return boxLabel + '-hingeLabel-' + i;
                });

            hingeLabel.enter()
                .append('text')
                .attr('class', 'hingeLabel boxplotValue');

            hingeLabel.transition()
                .duration(animationDuration)
                .attr('x', x(category) - hingeLength / 2)
                .attr('y', function (d) { return y(d); })
                .attr('dx', '-.35em')
                .attr('dy', '.35em')
                .attr('text-anchor', 'end')
                .text(function (d) { return d; });

            var connection = box.selectAll('.connection')
                .data([[params.upperWhisker, params.upperHinge], [params.lowerWhisker, params.lowerHinge]],
                    function (d, i) { return boxLabel + '-connection-' + i; });

            connection.enter()
                .append('line')
                .attr('class', 'connection');

            connection.transition()
                .duration(animationDuration)
                .attr('x1', x(category))
                .attr('y1', function (d) { return y(d[0]); })
                .attr('x2', x(category))
                .attr('y2', function (d) { return y(d[1]); });

            var upperBox = box.selectAll('.box.upper')
                .data(params.upperHinge, function (d) { return d; });

            upperBox.enter()
                .append('rect')
                .attr('class', 'box upper');

            upperBox.transition()
                .duration(animationDuration)
                .attr('x', x(category) - hingeLength / 2)
                .attr('y', y(params.upperHinge))
                .attr('height', Math.abs(y(params.upperHinge) - y(params.median)))
                .attr('width', hingeLength);

            var lowerBox = box.selectAll('.box.lower')
                .data(params.lowerHinge, function (d) { return d; });

            lowerBox.enter()
                .append('rect')
                .attr('class', 'box lower');

            lowerBox.transition()
                .duration(animationDuration)
                .attr('x', x(category) - hingeLength / 2)
                .attr('y', y(params.median))
                .attr('height', Math.abs(y(params.median) - y(params.lowerHinge)))
                .attr('width', hingeLength);

            var medianLabel = box.selectAll('.medianLabel')
                .data(params.median, function (d) { return d;});

            medianLabel.enter()
                .append('text')
                .attr('class', 'medianLabel boxplotValue');

            medianLabel.transition()
                .duration(animationDuration)
                .attr('x', x(category) + hingeLength / 2)
                .attr('y', function () { return y(params.median); })
                .attr('dx', '.35em')
                .attr('dy', '.35em')
                .attr('text-anchor', 'start')
                .text(params.median);

            var point = box.selectAll('.point')
                .data(params.rawData, function (d) { return boxLabel + '-' + d.patientID; });

            point.enter()
                .append('circle')
                .attr('cx', function (d) {
                    return jitterCheck.checked ? x(category) + boxWidth * jitterWidth * d.jitter : x(category);
                })
                .attr('r', 0)
                .attr('fill', function (d) { return colorScale(d.value); })
                .on('mouseover', function (d) {
                    tip.show('Value: ' + d.value + '</br>' +
                        'PatientID: ' + d.patientID + '</br>' +
                        'Outlier: ' + d.outlier);
                })
                .on('mouseout', function () {
                    tip.hide();
                });

            point
                .attr('class', function (d) {
                    return 'point patientID-' + smartRUtils.makeSafeForCSS(d.patientID) +
                        (d.outlier ? ' outlier ' : ' ') + smartRUtils.makeSafeForCSS(boxLabel);
                }) // This is here and not in the .enter() because points might become outlier on removal of other points
                .transition()
                .duration(animationDuration)
                .attr('cy', function (d) { return y(d.value); })
                .attr('r', 3);

            point.exit()
                .transition()
                .duration(animationDuration)
                .attr('r', 0)
                .remove();

            var yCopy = y.copy();
            yCopy.domain([params.lowerWhisker, params.upperWhisker]);
            var kde = kernelDensityEstimator(epanechnikovKernel(6), yCopy.ticks(100));
            var values = params.rawData.map(function (d) { return d.value; });
            var estFun = kde(values);
            var kdeDomain = d3.extent(estFun, function (d) { return d[1]; });
            var kdeScale = d3.scale.linear()
                .domain(kdeDomain)
                .range([0, boxWidth / 2]);
            var lineGen = d3.svg.line()
                .x(function (d) { return x(category) - kdeScale(d[1]); })
                .y(function (d) { return y(d[0]); });

            box.append('path')
                .attr('class', 'line')
                .attr('visibility', 'hidden')
                .datum(estFun)
                .transition()
                .duration(animationDuration)
                .attr('d', lineGen);
        }

        function removeBrush() {
            d3.selectAll('.brush')
                .call(brush.clear());
        }

        function reset() {
            removeBrush();
            excludedPatientIDs = [];
            currentSelection = [];
            excludeSelection(); // Abusing the method because I can
        }

    }

}]);

