//# sourceURL=d3Heatmap.js

'use strict';

window.smartRApp.directive('heatmapPlot', [
    'smartRUtils',
    '$rootScope',
    function(smartRUtils, $rootScope) {

        return {
            restrict: 'E',
            scope: {
                data: '=',
                width: '@',
                height: '@',
                params: '='
            },
            templateUrl: $rootScope.smartRPath +  '/js/smartR/_angular/templates/heatmap.html',
            link: function(scope, element) {
                var viz = element.children()[1];
                /**
                 * Watch data model (which is only changed by ajax calls when we want to (re)draw everything)
                 */
                scope.$watch('data', function(newValue) {
                    scope.showControls = false;
                    angular.element(viz).empty();
                    if (angular.isArray(newValue.fields)) {
                        scope.showControls = true;
                        createHeatmap(scope, viz);
                    }
                }, true);
            }
        };

        function createHeatmap(scope, root) {
            console.log(scope.data);
            var ANIMATION_DURATION = 1500;

            var fields = scope.data.fields;
            var extraFields = scope.data.extraFields;
            var features = scope.data.features.constructor === Array ? scope.data.features : [];

            var patientIDs = scope.data.patientIDs.map(function(d) { return parseInt(d); });
            var colNames = scope.data.colNames; // unique
            var rowNames = scope.data.rowNames; // unique

            var longestColName = colNames.reduce(function(prev, curr) {
                return curr.length > prev.length ? curr : prev;
            }, '');
            var longestColNameLength = smartRUtils.getTextWidth(longestColName);

            var longestRowName = rowNames.reduce(function(prev, curr) {
                return curr.length > prev.length ? curr : prev;
            }, '');
            var longestRowNameLength = smartRUtils.getTextWidth(longestRowName);

            var originalColNames = colNames.slice();
            var originalRowNames = rowNames.slice();

            var ranking = scope.data.ranking[0].toUpperCase();
            var statistics = scope.data.allStatValues;

            var maxRows = scope.data.maxRows[0];

            var geneCardsAllowed = JSON.parse(scope.params.geneCardsAllowed);

            var gridFieldWidth = 20;
            var gridFieldHeight = 10;
            var dendrogramHeight = 300;
            var histogramHeight = 200;
            var legendWidth = 200;
            var legendHeight = 40;

            var margin = {
                top: gridFieldHeight * 2 + features.length * gridFieldHeight + dendrogramHeight + 100,
                right: gridFieldWidth + 300 + dendrogramHeight,
                bottom: 10,
                left: histogramHeight
            };

            var width = gridFieldWidth * colNames.length;
            var height = gridFieldHeight * rowNames.length;

            // FIXME: This is here because the sizing of the whole heatmap is kind of messed up
            // At one point in the future we need to fix this
            smartRUtils.prepareWindowSize(width * 2 + margin.left + margin.right, height * 2 + margin.top + margin.right);

            var selectedColNames = [];

            var scale = null;
            var histogramScale = null;

            var animationCheck = smartRUtils.getElementWithoutEventListeners('sr-heatmap-animate-check');
            animationCheck.checked = fields.length < 10000;

            var zoomRange = smartRUtils.getElementWithoutEventListeners('sr-heatmap-zoom-range');
            zoomRange.addEventListener('mouseup', function() { zoom(parseInt(zoomRange.value)); });
            zoomRange.value = 100;

            var cutoffBtn = smartRUtils.getElementWithoutEventListeners('sr-heatmap-cutoff-btn');
            cutoffBtn.addEventListener('click', cutoff);

            var cutoffRange = smartRUtils.getElementWithoutEventListeners('sr-heatmap-cutoff-range');
            cutoffRange.addEventListener('mouseup', function() { animateCutoff(parseInt(cutoffRange.value)); });
            cutoffRange.setAttribute('max', maxRows);
            cutoffRange.value = 0;

            var clusterSelect = smartRUtils.getElementWithoutEventListeners('sr-heatmap-cluster-select');
            clusterSelect.addEventListener('change', function() { cluster(clusterSelect.value); });
            clusterSelect.disabled = maxRows < 2;
            clusterSelect.selectedIndex = 0;

            var clusterRowCheck = smartRUtils.getElementWithoutEventListeners('sr-heatmap-row-check');
            clusterRowCheck.disabled = maxRows < 2;
            clusterRowCheck.checked = true;

            var clusterColCheck = smartRUtils.getElementWithoutEventListeners('sr-heatmap-col-check');
            clusterColCheck.disabled = maxRows < 2;
            clusterColCheck.checked = true;

            var colorSelect = smartRUtils.getElementWithoutEventListeners('sr-heatmap-color-select');
            colorSelect.addEventListener('change', function() { updateColors(colorSelect.value); });
            colorSelect.selectedIndex = 0;

            var rankingSelect = smartRUtils.getElementWithoutEventListeners('sr-heatmap-ranking-select');
            rankingSelect.addEventListener('change', function() { setRanking(rankingSelect.value); });
            while (rankingSelect.firstChild) {
                rankingSelect.removeChild(rankingSelect.firstChild);
            }
            for (var stat in statistics[0]) { // collect existing statistics headers
                if (statistics[0].hasOwnProperty(stat) && stat !== 'ROWNAME') {
                    var option = document.createElement('option');
                    if (ranking === stat) {
                        option.selected = true;
                    }
                    option.setAttribute('value', stat);
                    option.innerHTML = stat.toLowerCase();
                    rankingSelect.appendChild(option);
                }
            }

            function setScales() {
                scale = d3.scale.linear()
                    .domain(d3.extent(statistics.map(function(d) { return d[ranking]; })))
                    .range((ranking === 'PVAL' || ranking === 'ADJPVAL') ? [histogramHeight, 0] : [0, histogramHeight]);

                histogramScale = function(value) {
                    return (ranking === 'TTEST' || ranking === 'LOGFOLD') ? scale(Math.abs(value)) : scale(value);
                };
            }
            setScales();

            function getInternalSortValue(value) {
                switch (ranking) {
                    case 'PVAL':
                    case 'ADJPVAL':
                        return 1 - value;
                    default:
                        return value;
                }
            }

            var heatmap = d3.select(root).append('svg')
                .attr('width', (width + margin.left + margin.right) * 4)
                .attr('height', (height + margin.top + margin.bottom) * 4)
                .attr('class', 'visualization')
                .append('g')
                .attr('transform', 'translate(' + margin.left + ',' + margin.top + ')');

            function adjustDimensions() {
                // gridFieldWidth/gridFieldHeight are adjusted outside as the zoom changes
                $(heatmap[0]).closest('svg')
                    .attr('width', margin.left + margin.right + (gridFieldWidth * colNames.length))
                    .attr('height', margin.top + margin.bottom + (gridFieldHeight * rowNames.length));
            }

            adjustDimensions();

            var tip = d3.tip()
                .attr('class', 'd3-tip')
                .offset([-10, 0])
                .html(function(d) { return d; });

            heatmap.call(tip);

            var featureItems = heatmap.append('g');
            var squareItems = heatmap.append('g');
            var colSortItems = heatmap.append('g');
            var selectItems = heatmap.append('g');
            var colNameItems = heatmap.append('g');
            var rowSortItems = heatmap.append('g');
            var significanceSortItems = heatmap.append('g');
            var labelItems = heatmap.append('g');
            var barItems = heatmap.append('g');
            var legendItems = heatmap.append('g');

            // this code is needed for the legend generation
            var zScores = fields.map(function(d) { return d.ZSCORE; });
            var maxZScore = zScores.reduce(function(prev, curr) { return curr > prev ? curr : prev; });
            var minZScore = zScores.reduce(function(prev, curr) { return curr < prev ? curr : prev; });
            var steps = [];
            for (var i = minZScore; i < maxZScore; i+= (maxZScore - minZScore) / 50) {
                steps.push(i);
            }

            function updateHeatmap() {
                updateHeatmapTable();
                var square = squareItems.selectAll('.square')
                    .data(fields);

                square.enter()
                    .append('rect')
                    .attr('class', function(d) {
                        return 'square colname-' + smartRUtils.makeSafeForCSS(d.COLNAME) +
                            ' rowname-' + smartRUtils.makeSafeForCSS(d.ROWNAME);
                    })
                    .attr('x', function(d) { return colNames.indexOf(d.COLNAME) * gridFieldWidth; })
                    .attr('y', function(d) { return rowNames.indexOf(d.ROWNAME) * gridFieldHeight; })
                    .attr('width', gridFieldWidth)
                    .attr('height', gridFieldHeight)
                    .on('mouseover', function(d) {
                        d3.select('.colname.colname-' + smartRUtils.makeSafeForCSS(d.COLNAME))
                            .classed('highlight', true);
                        d3.select('.rowname.rowname-' + smartRUtils.makeSafeForCSS(d.ROWNAME))
                            .classed('highlight', true);

                        var html = 'Log2: ' + d.VALUE + '<br/>' +
                            'z-Score: ' + d.ZSCORE + '<br/>' +
                            'Column: ' + d.COLNAME + '<br/>' +
                            'Row: ' + d.ROWNAME + '<br/>' +
                            'PatientId: ' + d.PATIENTID + '</br>' +
                            'Subset: ' + d.SUBSET + '<br/>';

                        tip.show(html);
                    })
                    .on('mouseout', function() {
                        d3.selectAll('.colname').classed('highlight', false);
                        d3.selectAll('.rowname').classed('highlight', false);
                        tip.hide();
                    });

                square.transition()
                    .duration(animationCheck.checked ? ANIMATION_DURATION : 0)
                    .attr('x', function(d) { return colNames.indexOf(d.COLNAME) * gridFieldWidth; })
                    .attr('y', function(d) { return rowNames.indexOf(d.ROWNAME) * gridFieldHeight; })
                    .attr('width', gridFieldWidth)
                    .attr('height', gridFieldHeight);

                var colSortText = colSortItems.selectAll('.colSortText')
                    .data(colNames);

                colSortText.enter()
                    .append('text')
                    .attr('class', 'colSortText')
                    .attr('dy', '0.35em')
                    .style('text-anchor', 'middle')
                    .text('↑↓');

                colSortText.transition()
                    .duration(animationCheck.checked ? ANIMATION_DURATION : 0)
                    .attr('x', function(d, i) { return i * gridFieldWidth + 0.5 * gridFieldWidth; })
                    .attr('y', -2 - gridFieldHeight + 0.5 * gridFieldHeight)
                    .style('font-size', gridFieldHeight + 'px');

                function getValueForSquareSorting(colName, rowName) {
                    var square = d3.select('.square' + '.colname-' + smartRUtils.makeSafeForCSS(colName) +
                        '.rowname-' + smartRUtils.makeSafeForCSS(rowName));
                    return square[0][0] ? square.property('__data__').ZSCORE : Number.NEGATIVE_INFINITY;
                }

                function isSorted(arr) {
                    return arr.every(function(d, i) {
                        return i === arr.length - 1 || arr[i][1] >= arr[i + 1][1];
                    });
                }

                var colSortBox = colSortItems.selectAll('.colSortBox')
                    .data(colNames);

                colSortBox.enter()
                    .append('rect')
                    .attr('class', function(d, i) { return 'box colSortBox idx-' + i; })
                    .on('click', function(colName) {
                        var rowValues = rowNames.map(function(rowName, idx) {
                            return [idx, getValueForSquareSorting(colName, rowName)];
                        });
                        if (isSorted(rowValues)) {
                            rowValues.sort(function(a, b) { return a[1] - b[1]; });
                        } else {
                            rowValues.sort(function(a, b) { return b[1] - a[1]; });
                        }
                        var sortValues = rowValues.map(function(rowValue) { return rowValue[0]; });
                        updateRowOrder(sortValues);

                        d3.selectAll('.colSortBox').classed('sortedBy', false);
                        d3.selectAll('.significanceSortBox').classed('sortedBy', false);
                        d3.select(this).classed('sortedBy', true);

                        var sortedByBoxes = d3.selectAll('.sortedBy').filter('.rowSortBox');
                        if (!sortedByBoxes.empty()) {
                            var className = sortedByBoxes.attr('class');
                            var sortedByIdx = parseInt(className.match(/idx-(\d+)/)[1]);
                            d3.selectAll('.rowSortBox').classed('sortedBy', false);
                            d3.selectAll('.rowSortBox').filter('.idx-' + sortValues.indexOf(sortedByIdx)).classed('sortedBy', true);
                        }
                    });

                colSortBox.transition()
                    .duration(animationCheck.checked ? ANIMATION_DURATION : 0)
                    .attr('x', function(d, i) { return i * gridFieldWidth; })
                    .attr('y', -2 - gridFieldHeight)
                    .attr('width', gridFieldWidth)
                    .attr('height', gridFieldHeight);

                var rowSortText = rowSortItems.selectAll('.rowSortText')
                    .data(rowNames);

                rowSortText.enter()
                    .append('text')
                    .attr('class', 'rowSortText')
                    .attr('dy', '0.35em')
                    .style('text-anchor', 'middle')
                    .text('↑↓');

                rowSortText.transition()
                    .duration(animationCheck.checked ? ANIMATION_DURATION : 0)
                    .style('font-size', gridFieldHeight + 'px')
                    .attr('transform', function(d, i) {
                        return 'translate(' + (width + 2 + 0.5 * gridFieldWidth) + ',0)' + 'translate(0,' +
                            (i * gridFieldHeight + 0.5 * gridFieldHeight) + ')rotate(-90)';
                    });

                var rowSortBox = rowSortItems.selectAll('.rowSortBox')
                    .data(rowNames);

                rowSortBox.enter()
                    .append('rect')
                    .attr('class', function(d, i) { return 'box rowSortBox idx-' + i; })
                    .on('click', function(rowName) {
                        var colValues = colNames.map(function(colName, idx) {
                            return [idx, getValueForSquareSorting(colName, rowName)];
                        });
                        if (isSorted(colValues)) {
                            colValues.sort(function(a, b) { return a[1] - b[1]; });
                        } else {
                            colValues.sort(function(a, b) { return b[1] - a[1]; });
                        }
                        var sortValues = colValues.map(function(colValue) { return colValue[0]; });
                        updateColOrder(sortValues);

                        d3.selectAll('.rowSortBox').classed('sortedBy', false);
                        d3.selectAll('.featureSortBox').classed('sortedBy', false);
                        d3.select(this).classed('sortedBy', true);
                        var sortedByBoxes = d3.selectAll('.sortedBy').filter('.colSortBox');
                        if (!sortedByBoxes.empty()) {
                            var className = sortedByBoxes.attr('class');
                            var sortedByIdx = parseInt(className.match(/idx-(\d+)/)[1]);
                            d3.selectAll('.colSortBox').classed('sortedBy', false);
                            d3.selectAll('.colSortBox').filter('.idx-' + sortValues.indexOf(sortedByIdx)).classed('sortedBy', true);
                        }
                    });

                rowSortBox.transition()
                    .duration(animationCheck.checked ? ANIMATION_DURATION : 0)
                    .attr('x', width + 2)
                    .attr('y', function(d, i) { return i * gridFieldHeight; })
                    .attr('width', gridFieldWidth)
                    .attr('height', gridFieldHeight);

                var significanceSortText = significanceSortItems.selectAll('.significanceSortText')
                    .data(['something']);

                significanceSortText.enter()
                    .append('text')
                    .attr('class', 'significanceSortText')
                    .attr('dy', '0.35em')
                    .style('text-anchor', 'middle')
                    .text('↑↓');

                significanceSortText.transition()
                    .duration(animationCheck.checked ? ANIMATION_DURATION : 0)
                    .style('font-size', gridFieldHeight + 'px')
                    .attr('x', -gridFieldWidth - 10 + 0.5 * gridFieldWidth)
                    .attr('y', -2 - gridFieldHeight + 0.5 * gridFieldHeight);

                var significanceSortBox = significanceSortItems.selectAll('.significanceSortBox')
                    .data(['something']);

                significanceSortBox.enter()
                    .append('rect')
                    .attr('class', 'box significanceSortBox')
                    .on('click', function() {
                        var rowValues = statistics.map(function(d) { return d[ranking]; })
                            .map(function(significanceValue, idx) {
                            return [idx, getInternalSortValue(significanceValue)];
                        });

                        if (isSorted(rowValues)) {
                            rowValues.sort(function(a, b) { return a[1] - b[1]; });
                        } else {
                            rowValues.sort(function(a, b) { return b[1] - a[1]; });
                        }
                        var sortValues = rowValues.map(function(rowValue) { return rowValue[0]; });
                        updateRowOrder(sortValues);

                        d3.selectAll('.colSortBox').classed('sortedBy', false);
                        d3.select(this).classed('sortedBy', true);

                        var sortedByBoxes = d3.selectAll('.sortedBy').filter('.rowSortBox');
                        if (!sortedByBoxes.empty()) {
                            var className = sortedByBoxes.attr('class');
                            var sortedByIdx = parseInt(className.match(/idx-(\d+)/)[1]);
                            d3.selectAll('.rowSortBox').classed('sortedBy', false);
                            d3.selectAll('.rowSortBox').filter('.idx-' + sortValues.indexOf(sortedByIdx)).classed('sortedBy', true);
                        }
                    });

                significanceSortBox.transition()
                    .duration(animationCheck.checked ? ANIMATION_DURATION : 0)
                    .attr('x', -gridFieldWidth - 10)
                    .attr('y', -2 - gridFieldHeight)
                    .attr('width', gridFieldWidth)
                    .attr('height', gridFieldHeight);

                var selectText = selectItems.selectAll('.selectText')
                    .data(colNames, function(d) {
                        return d;
                    });

                selectText.enter()
                    .append('text')
                    .attr('class', function(d) { return 'selectText colname-' + smartRUtils.makeSafeForCSS(d); })
                    .attr('dy', '0.35em')
                    .style('text-anchor', 'middle')
                    .text('□');

                selectText.transition()
                    .duration(animationCheck.checked ? ANIMATION_DURATION : 0)
                    .style('font-size', gridFieldHeight + 'px')
                    .attr('x', function(d, i) { return i * gridFieldWidth + 0.5 * gridFieldWidth; })
                    .attr('y', -2 - gridFieldHeight * 2 + 0.5 * gridFieldHeight);

                var selectBox = selectItems.selectAll('.selectBox')
                    .data(colNames);

                selectBox.enter()
                    .append('rect')
                    .attr('class', 'box selectBox')
                    .on('click', function(colName) { selectCol(colName); });

                selectBox.transition()
                    .duration(animationCheck.checked ? ANIMATION_DURATION : 0)
                    .attr('x', function(d, i) { return i * gridFieldWidth; })
                    .attr('y', -2 - gridFieldHeight * 2)
                    .attr('width', gridFieldWidth)
                    .attr('height', gridFieldHeight);

                var colName = colNameItems.selectAll('.colname')
                    .data(colNames.map(function(colName, i) {
                        return { colName: colName, patientID: patientIDs[i] };
                    }), function(d) { return d.colName; });

                colName.enter()
                    .append('text')
                    .attr('class', function(d) {
                        return 'colname colname-' + smartRUtils.makeSafeForCSS(d.colName) +
                            ' patientID-' + smartRUtils.makeSafeForCSS(d.patientID);
                    })
                    .style('text-anchor', 'start')
                    .text(function(d) { return d.colName; });

                colName.transition()
                    .duration(animationCheck.checked ? ANIMATION_DURATION : 0)
                    .style('font-size', gridFieldHeight + 'px')
                    .attr('transform', function(d) {
                        return 'translate(' + (colNames.indexOf(d.colName) * gridFieldWidth) + ',0)' +
                            'translate(' + (gridFieldWidth / 2) + ',' + (-4 - gridFieldHeight * 2) + ')rotate(-45)';
                    });

                var rowName = labelItems.selectAll('.rowname')
                    .data(rowNames, function(d) { return d; });

                rowName.enter()
                    .append('text')
                    .attr('class', function(d) { return 'rowname rowname-' + smartRUtils.makeSafeForCSS(d); })
                    .attr('dy', '0.35em')
                    .style('text-anchor', 'start')
                    .text(function(d) { return d; })
                    .on('click', function(d) {
                        var genes = d.split('--');
                        genes.shift();
                        var urls = [];
                        if (geneCardsAllowed) {
                            genes.forEach(function(gene) {
                                urls.push('http://www.genecards.org/cgi-bin/carddisp.pl?gene=' + gene);
                            });
                        } else {
                            genes.forEach(function(gene) {
                                urls.push('https://www.ebi.ac.uk/ebisearch/search.ebi?db=allebi&query=' + gene);
                            });
                        }
                        urls.forEach(function(url) {
                            window.open(url);
                        });
                    });

                rowName.transition()
                    .duration(animationCheck.checked ? ANIMATION_DURATION : 0)
                    .style('font-size', gridFieldHeight + 'px')
                    .attr('x', width + gridFieldWidth + 7)
                    .attr('y', function(d) { return rowNames.indexOf(d) * gridFieldHeight + 0.5 * gridFieldHeight; });

                var bar = barItems.selectAll('.bar')
                    .data(statistics, function(d, i) { return i; });

                bar.enter()
                    .append('rect')
                    .on('mouseover', function(d) {
                        var html = '';
                        for (var key in d) {
                            if (d.hasOwnProperty(key) && key !== 'ROWNAME') {
                                html += (key === ranking ? '(ranked by) ' : '') + key + ': ' + d[key] + '<br/>';
                            }
                        }
                        tip.show(html);
                        d3.selectAll('.square.rowname-' + smartRUtils.makeSafeForCSS(d.ROWNAME))
                            .classed('squareHighlighted', true);
                        d3.select('.rowname.rowname-' + smartRUtils.makeSafeForCSS(d.ROWNAME))
                            .classed('highlight', true);
                    })
                    .on('mouseout', function() {
                        tip.hide();
                        d3.selectAll('.square').classed('squareHighlighted', false);
                        d3.selectAll('.rowname').classed('highlight', false);
                    });

                bar.attr('class', function(d) { return 'bar rowname-' + smartRUtils.makeSafeForCSS(d.ROWNAME); })
                    .transition()
                    .duration(animationCheck.checked ? ANIMATION_DURATION : 0)
                    .attr('width', function(d) { return histogramScale(d[ranking]); })
                    .attr('height', gridFieldHeight)
                    .attr('x', function(d) { return -histogramScale(d[ranking]); })
                    .attr('y', function(d) { return gridFieldHeight * rowNames.indexOf(d.ROWNAME); })
                    .style('fill', function(d) { return d[ranking] > 0 ? '#990000' : 'steelblue'; });

                var featurePosY = -gridFieldWidth * 2 - longestColNameLength + 20;

                var extraSquare = featureItems.selectAll('.extraSquare')
                    .data(extraFields);

                extraSquare.enter()
                    .append('rect')
                    .attr('class', function(d) {
                        return 'extraSquare patientID-' + smartRUtils.makeSafeForCSS(d.PATIENTID) +
                            ' rowname-' + smartRUtils.makeSafeForCSS(d.ROWNAME);
                    })
                    .on('mouseover', function(d) {
                        d3.select('.colname.patientID-' + smartRUtils.makeSafeForCSS(d.PATIENTID)).classed('highlight', true);
                        d3.select('.feature.feature-' + smartRUtils.makeSafeForCSS(d.ROWNAME)).classed('highlight', true);
                        var html = 'Value: ' + d.VALUE + '<br/>' +
                            (d.ZSCORE ? 'z-Score: ' + d.ZSCORE + '<br/>' : '') +
                            'Column: ' + d.COLNAME + '<br/>' +
                            'Row: ' + d.ROWNAME + '<br/>' +
                            'PatientId: ' + d.PATIENTID + '<br/>' +
                            'Type: ' + d.TYPE + '<br/>' +
                            'Subset: ' + d.SUBSET;

                        tip.show(html);
                    })
                    .on('mouseout', function() {
                        tip.hide();
                        d3.selectAll('.colname').classed('highlight', false);
                        d3.selectAll('.feature').classed('highlight', false);
                    });

                extraSquare.transition()
                    .duration(animationCheck.checked ? ANIMATION_DURATION : 0)
                    .attr('x', function(d) { return patientIDs.indexOf(d.PATIENTID) * gridFieldWidth; })
                    .attr('y', function(d) { return featurePosY - features.indexOf(d.ROWNAME) * gridFieldHeight; })
                    .attr('width', gridFieldWidth)
                    .attr('height', gridFieldHeight);

                var feature = featureItems.selectAll('.feature')
                    .data(features);

                feature.enter()
                    .append('text')
                    .attr('class', function(d) { return 'feature feature-' + smartRUtils.makeSafeForCSS(d); })
                    .attr('dy', '0.35em')
                    .style('text-anchor', 'start')
                    .text(function(d) {
                        return d === 'Cohort' ? d : smartRUtils.shortenConcept(d);
                    });

                feature.transition()
                    .duration(animationCheck.checked ? ANIMATION_DURATION : 0)
                    .style('font-size', gridFieldHeight + 'px')
                    .attr('x', width + gridFieldWidth + 7)
                    .attr('y', function(d) { return featurePosY - features.indexOf(d) * gridFieldHeight + gridFieldHeight / 2; });

                var featureSortText = featureItems.selectAll('.featureSortText')
                    .data(features);

                featureSortText.enter()
                    .append('text')
                    .attr('class', 'featureSortText')
                    .attr('dy', '0.35em')
                    .style('text-anchor', 'middle')
                    .text('↑↓');

                featureSortText.transition()
                    .duration(animationCheck.checked ? ANIMATION_DURATION : 0)
                    .style('font-size', gridFieldHeight + 'px')
                    .attr('transform', function(d) {
                        return 'translate(' + (width + 2 + 0.5 * gridFieldWidth) + ',0)' + 'translate(0,' +
                            (featurePosY - features.indexOf(d) * gridFieldHeight + gridFieldHeight / 2) + ')rotate(-90)';
                    });

                var featureSortBox = featureItems.selectAll('.featureSortBox')
                    .data(features);

                featureSortBox.enter()
                    .append('rect')
                    .attr('class', 'box featureSortBox')
                    .on('click', function(feature) {
                        var featureValues = patientIDs.map(function(patientID, idx) {
                            var css = 'extraSquare rowname-' + smartRUtils.makeSafeForCSS(feature) +
                                ' patientID-' + smartRUtils.makeSafeForCSS(patientID);
                            var elements = document.getElementsByClassName(css);
                            var value = (-Math.pow(2, 32)).toString(); // if square does not exist
                            if (elements.length > 0) {
                                var square = d3.select(elements[0]);
                                value = square.property('__data__').VALUE;
                            }
                            return [idx, value];
                        });
                        if (isSorted(featureValues)) {
                            featureValues.sort(function(a, b) {
                                var diff = a[1] - b[1];
                                return isNaN(diff) ? a[1].localeCompare(b[1]) : diff;
                            });
                        } else {
                            featureValues.sort(function(a, b) {
                                var diff = b[1] - a[1];
                                return isNaN(diff) ? b[1].localeCompare(a[1]) : diff;
                            });
                        }
                        var sortValues = featureValues.map(function(featureValue) {
                            return featureValue[0];
                        });
                        updateColOrder(sortValues);

                        d3.selectAll('.rowSortBox').classed('sortedBy', false);
                        d3.selectAll('.featureSortBox').classed('sortedBy', false);
                        d3.select(this).classed('sortedBy', true);
                        var sortedByBoxes = d3.selectAll('.sortedBy').filter('.colSortBox');
                        if (!sortedByBoxes.empty()) {
                            var className = sortedByBoxes.attr('class');
                            var sortedByIdx = parseInt(className.match(/idx-(\d+)/)[1]);
                            d3.selectAll('.colSortBox').classed('sortedBy', false);
                            d3.selectAll('.colSortBox').filter('.idx-' + sortValues.indexOf(sortedByIdx)).classed('sortedBy', true);
                        }
                    });


                featureSortBox.transition()
                    .duration(animationCheck.checked ? ANIMATION_DURATION : 0)
                    .attr('x', width + 2)
                    .attr('y', function(d) { return featurePosY - features.indexOf(d) * gridFieldHeight; })
                    .attr('width', gridFieldWidth)
                    .attr('height', gridFieldHeight);
            }

            function updateHeatmapTable() {
                d3.selectAll('.sr-heatmap-table').remove();

                var HEADER = ['ROWNAME'];
                for (var stat in statistics[0]) { // collect existing statistics headers
                    if (statistics[0].hasOwnProperty(stat) && stat !== 'ROWNAME') {
                        HEADER.push(stat);
                    }
                }
                var table = d3.select(root).append('table')
                    .attr('class', 'sr-heatmap-table');
                var thead = table.append('thead');
                var tbody = table.append('tbody');

                thead.append('tr')
                    .selectAll('th')
                    .data(HEADER)
                    .enter()
                    .append('th')
                    .text(function(d) {
                        return d;
                    });

                var probeIDs = [];
                var entities = [];
                rowNames.forEach(function(rowName) {
                    probeIDs.push(rowName.match(/.+(?=--)/)[0]);
                    entities.push(rowName.match(/.+?--(.*)/)[1]);
                });

                var rows = tbody.selectAll('tr')
                    .data(statistics)
                    .enter()
                    .append('tr');

                rows.selectAll('td')
                    .data(function(d, i) {
                        return HEADER.map(function(column) {
                            return {column: column, value: statistics[i][column]};
                        });
                    })
                    .enter()
                    .append('td')
                    .text(function(d) {
                        return d.value;
                    });
            }

            function zoom(zoomLevel) {
                zoomLevel /= 100;
                gridFieldWidth = 20 * zoomLevel;
                gridFieldHeight = 10 * zoomLevel;
                width = gridFieldWidth * colNames.length;
                height = gridFieldHeight * rowNames.length;
                heatmap
                    .attr('width', width + margin.left + margin.right)
                    .attr('height', width + margin.top + margin.bottom);
                var temp = ANIMATION_DURATION;
                ANIMATION_DURATION = 0;
                updateHeatmap();
                reloadDendrograms();
                ANIMATION_DURATION = temp;
                adjustDimensions();
            }

            var cutoffLevel = 0;

            function animateCutoff(cutoff) {
                cutoff = Math.floor(cutoff);
                cutoffLevel = cutoff;
                d3.selectAll('.square')
                    .classed('cuttoffHighlight', false);
                d3.selectAll('.bar')
                    .classed('cuttoffHighlight', false);
                statistics.map(function(d) { return d[ranking]; })
                    .sort(function(a, b) { return a - b; })
                    .filter(function(d, i) { return i < cutoff; })
                    .forEach(function(d) {
                        d3.select('.bar.idx-' + smartRUtils.makeSafeForCSS(d[0])).classed('cuttoffHighlight', true);
                        d3.selectAll('.square.rowname-' + smartRUtils.makeSafeForCSS(rowNames[d[0]])).classed('cuttoffHighlight', true);
                    });
            }

            function cutoff() {
                //HeatmapService.startScriptExecution({
                //    taskType: 'run',
                //    arguments: params,
                //    onUltimateSuccess: HeatmapService.runAnalysisSuccess,
                //    onUltimateFailure: HeatmapService.runAnalysisFailed,
                //    phase: 'run',
                //    progressMessage: 'Calculating',
                //    successMessage: undefined
                //});
                // TODO: Use ajax service to be provided by ajaxServices.js to re-compute analysis
                // with new arguments (in this case filter for cut-off)
                scope.params.max_row = maxRows - cutoffLevel - 1;
                $('run-button input').click();
            }

            function reloadDendrograms() {
                if (colDendrogramVisible) {
                    removeColDendrogram();
                    createColDendrogram();
                }
                if (rowDendrogramVisible) {
                    removeRowDendrogram();
                    createRowDendrogram();
                }
            }

            function selectCol(colName) {
                var colSquares = d3.selectAll('.square.colname-' + smartRUtils.makeSafeForCSS(colName));
                if (colSquares.classed('selected')) {
                    var index = selectedColNames.indexOf(colName);
                    selectedColNames.splice(index, 1);
                    colSquares
                        .classed('selected', false);
                    d3.select('.selectText.colname-' + smartRUtils.makeSafeForCSS(colName))
                        .text('□');
                } else {
                    selectedColNames.push(colName);
                    colSquares.classed('selected', true);
                    d3.select('.selectText.colname-' + smartRUtils.makeSafeForCSS(colName))
                        .text('■');
                }
                if (selectedColNames.length !== 0) {
                    d3.selectAll('.square:not(.selected)')
                        .attr('opacity', 0.4);
                } else {
                    d3.selectAll('.square')
                        .attr('opacity', 1);
                }
            }

            var colorScale;
            function updateColors(schema) {
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

                var redBlueScale = d3.scale.quantile()
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
                            colorSet.push(d3.rgb(0, 0, (255 * (NUM - i)) / NUM));
                        }
                        return colorSet.reverse();
                    }());

                var blueScale = d3.scale.linear()
                    .domain([0, 1])
                    .range(['#0000ff', '#e5e5ff']);

                var greenScale = d3.scale.linear()
                    .domain([0, 1])
                    .range(['#00ff00', '#e5ffe5']);

                var colorSchemas = {
                    redGreen: redGreenScale,
                    blueScale: blueScale,
                    redBlue: redBlueScale,
                    greenScale: greenScale
                };

                colorScale = colorSchemas[schema];

                d3.selectAll('.square')
                    .transition()
                    .duration(animationCheck.checked ? ANIMATION_DURATION : 0)
                    .attr('fill', function(d) {
                        return colorScale(1 / (1 + Math.pow(Math.E, -d.ZSCORE)));
                    });
                    
                d3.selectAll('.legendColor')
                    .transition()
                    .duration(animationCheck.checked ? ANIMATION_DURATION : 0)
                    .attr('fill', function(d) {
                        return colorScale(1 / (1 + Math.pow(Math.E, -d)));
                    });

                var featureColorSetBinary = ['#FF8000', '#FFFF00'];
                var featureColorSetSequential = [
                    'rgb(247,252,253)', 'rgb(224,236,244)', 'rgb(191,211,230)',
                    'rgb(158,188,218)', 'rgb(140,150,198)', 'rgb(140,107,177)',
                    'rgb(136,65,157)', 'rgb(129,15,124)', 'rgb(77,0,75)'
                ];
                var featureColorCategorical = d3.scale.category10();

                features.forEach(function(feature) {
                    d3.selectAll('.extraSquare.rowname-' + smartRUtils.makeSafeForCSS(feature))
                        .style('fill', function(d) {
                            switch (d.TYPE) {
                                case 'binary':
                                    return featureColorSetBinary[d.VALUE];
                                case 'subset':
                                    return featureColorSetBinary[d.VALUE - 1];
                                case 'numeric':
                                    colorScale.range(featureColorSetSequential);
                                    return colorScale(1 / (1 + Math.pow(Math.E, -d.ZSCORE)));
                                default:
                                    return featureColorCategorical(d.VALUE);
                            }
                        });
                });

                updateLegend();
            }

            function updateLegend() {
                var legendElementWidth = legendWidth / steps.length;
                var legendElementHeight = legendHeight;

                var legendColor = legendItems.selectAll('.legendColor')
                    .data(steps, function(d) { return d; });

                legendColor.enter()
                    .append('rect')
                    .attr('class', 'legendColor')
                    .attr('x', function(d, i) { return 5 - margin.left + i * legendElementWidth; })
                    .attr('y', 8 - margin.top + 100)
                    .attr('width', Math.ceil(legendElementWidth))
                    .attr('height', legendElementHeight);

                d3.selectAll('.legendText').remove();

                legendItems.append('text')
                    .attr('class', 'legendText')
                    .attr('x', 5 - margin.left)
                    .attr('y', 8 - margin.top + 100)
                    .attr('text-anchor', 'start')
                    .text(Math.min.apply(null, steps).toFixed(1));

                legendItems.append('text')
                    .attr('class', 'legendText')
                    .attr('x', 5 - margin.left + legendWidth)
                    .attr('y', 8 - margin.top + 100)
                    .attr('text-anchor', 'end')
                    .text(Math.max.apply(null, steps).toFixed(1));
            }

            function unselectAll() {
                d3.selectAll('.selectText')
                    .text('□');
                d3.selectAll('.square')
                    .classed('selected', false)
                    .attr('opacity', 1);
                selectedColNames = [];
            }

            var colDendrogramVisible = false;
            var colDendrogram;

            function createColDendrogram() {
                var w = 200;
                var colDendrogramWidth = gridFieldWidth * colNames.length;
                var spacing = gridFieldWidth * 2 + longestColNameLength + features.length * gridFieldHeight - 20;

                var cluster = d3.layout.cluster()
                    .size([colDendrogramWidth, w])
                    .separation(function() { return 1; });

                var diagonal = d3.svg.diagonal()
                    .projection(function(d) { return [d.x, -spacing - w + d.y]; });

                var colDendrogramNodes = cluster.nodes(colDendrogram);
                var colDendrogramLinks = cluster.links(colDendrogramNodes);

                heatmap.selectAll('.colDendrogramLink')
                    .data(colDendrogramLinks)
                    .enter().append('path')
                    .attr('class', 'colDendrogram link')
                    .attr('d', diagonal);

                heatmap.selectAll('.colDendrogramNode')
                    .data(colDendrogramNodes)
                    .enter().append('circle')
                    .attr('class', 'colDendrogram node')
                    .attr('r', 4.5)
                    .attr('transform', function(d) {
                        return 'translate(' + d.x + ',' + (-spacing - w + d.y) + ')';
                    }).on('click', function(d) {
                        var previousSelection = selectedColNames.slice();
                        unselectAll();
                        var leafs = d.index.split(' ');
                        for (var i = 0; i < leafs.length; i++) {
                            var colName = colNames[leafs[i]];
                            selectCol(colName);
                        }
                        if (previousSelection.sort().toString() === selectedColNames.sort().toString()) {
                            unselectAll();
                        }
                    })
                    .on('mouseover', function(d) {
                        tip.show('Height: ' + d.height);
                    })
                    .on('mouseout', function() {
                        tip.hide();
                    });
                colDendrogramVisible = true;
            }

            var rowDendrogramVisible = false;
            var rowDendrogram;

            function createRowDendrogram() {
                var h = 280;
                var rowDendrogramHeight = gridFieldHeight * rowNames.length;
                var spacing = gridFieldWidth + longestRowNameLength + 20;

                var cluster = d3.layout.cluster()
                    .size([rowDendrogramHeight, h])
                    .separation(function() {
                        return 1;
                    });

                var diagonal = d3.svg.diagonal()
                    .projection(function(d) {
                        return [width + spacing + h - d.y, d.x];
                    });

                var rowDendrogramNodes = cluster.nodes(rowDendrogram);
                var rowDendrogramLinks = cluster.links(rowDendrogramNodes);

                heatmap.selectAll('.rowDendrogramLink')
                    .data(rowDendrogramLinks)
                    .enter().append('path')
                    .attr('class', 'rowDendrogram link')
                    .attr('d', diagonal);

                heatmap.selectAll('.rowDendrogramNode')
                    .data(rowDendrogramNodes)
                    .enter().append('circle')
                    .attr('class', 'rowDendrogram node')
                    .attr('r', 4.5)
                    .attr('transform', function(d) {
                        return 'translate(' + (width + spacing + h - d.y) + ',' + d.x + ')';
                    }).on('click', function(d) {
                        var leafs = d.index.split(' ');
                        var genes = [];
                        leafs.forEach(function(leaf) {
                            var rowName = rowNames[leaf];
                            var split = rowName.split("--");
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
                    })
                    .on('mouseover', function(d) {
                        tip.show('Height: ' + d.height);
                    })
                    .on('mouseout', function() {
                        tip.hide();
                    });
                rowDendrogramVisible = true;
            }

            function removeColDendrogram() {
                heatmap.selectAll('.colDendrogram').remove();
                colDendrogramVisible = false;
            }

            function removeRowDendrogram() {
                heatmap.selectAll('.rowDendrogram').remove();
                rowDendrogramVisible = false;
            }

            function updateColOrder(sortValues, update) {
                update = typeof update === 'undefined' ? true : update;
                var newColnames = [];
                var newPatientIDs = [];
                sortValues.forEach(function(sortValue) {
                    newColnames.push(colNames[sortValue]);
                    newPatientIDs.push(patientIDs[sortValue]);
                });
                colNames = newColnames;
                patientIDs = newPatientIDs;
                unselectAll();
                removeColDendrogram();
                if (update) {
                    updateHeatmap();
                }
            }

            function updateRowOrder(sortValues, update) {
                update = typeof update === 'undefined' ? true : update;
                var sortedRowNames = [];
                var sortedStatistics = [];

                sortValues.forEach(function(sortValue) {
                    sortedRowNames.push(rowNames[sortValue]);
                    sortedStatistics.push(statistics[sortValue]);
                });
                rowNames = sortedRowNames;
                statistics = sortedStatistics;

                removeRowDendrogram();
                if (update) {
                    updateHeatmap();
                }
                animateCutoff();
            }

            function transformClusterOrderWRTInitialOrder(clusterOrder, initialOrder) {
                return clusterOrder.map(function(d) {
                    return initialOrder.indexOf(d);
                });
            }

            function getInitialRowOrder() {
                return rowNames.map(function(rowName) {
                    return originalRowNames.indexOf(rowName);
                });
            }

            function getInitialColOrder() {
                return colNames.map(function(colName) {
                    return originalColNames.indexOf(colName);
                });
            }

            var lastUsedClustering = null;

            function cluster(clustering) {
                if (!lastUsedClustering && typeof clustering === 'undefined') {
                    return; // Nothing should be done if clustering switches are turned on without clustering type set.
                }
                d3.selectAll('.box').classed('sortedBy', false);
                clustering = (typeof clustering === 'undefined') ? lastUsedClustering : clustering;
                var clusterData = scope.data.hclust[clustering];
                if (document.getElementById('sr-heatmap-row-check').checked && rowNames.length > 0) {
                    rowDendrogram = JSON.parse(clusterData[3]);
                    updateRowOrder(transformClusterOrderWRTInitialOrder(clusterData[1], getInitialRowOrder()), false);
                    createRowDendrogram(rowDendrogram);
                } else {
                    removeRowDendrogram();
                }
                if (document.getElementById('sr-heatmap-col-check').checked && colNames.length > 0) {
                    colDendrogram = JSON.parse(clusterData[2]);
                    updateColOrder(transformClusterOrderWRTInitialOrder(clusterData[0], getInitialColOrder()), false);
                    createColDendrogram(colDendrogram);
                } else {
                    removeColDendrogram();
                }
                lastUsedClustering = clustering;
                updateHeatmap();
            }

            function setRanking(method) {
                ranking = method;
                setScales();
                updateHeatmap();
            }

            function init() {
                var temp = ANIMATION_DURATION;
                ANIMATION_DURATION = 0;
                updateHeatmap();
                reloadDendrograms();
                updateColors('redGreen');
                updateColors('redGreen');
                ANIMATION_DURATION = temp;
            }

            init();

        }
    }]);
