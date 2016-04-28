//# sourceURL=d3Heatmap.js

'use strict';

window.smartRApp.directive('heatmapPlot', [
    'smartRUtils',
    'controlElements',
    function(smartRUtils, controlElements) {

        return {
            restrict: 'E',
            scope: {
                data: '=',
                width: '@',
                height: '@',
                params: '='
            },
            link: function(scope, element) {

                /**
                 * Watch data model (which is only changed by ajax calls when we want to (re)draw everything)
                 */
                scope.$watch('data', function(newValue) {
                    $(element[0]).empty();
                    if (angular.isArray(newValue.fields)) {
                        createHeatmap(scope.data, element[0], scope.params);
                    }
                }, true);
            }
        };

        function createHeatmap(data, root, params) {
            var animationDuration = 1500;
            var extraFields = data.extraFields === undefined ? [] : data.extraFields;
            var features = data.features === undefined ? [] : data.features;
            var fields = data.fields;
            var significanceValues = data.significanceValues;
            var logfoldValues = data.logfoldValues;
            var ttestValues = data.ttestValues;
            var pvalValues = data.pvalValues;
            var adjpvalValues = data.adjpvalValues;
            var bvalValues = data.bvalValues;
            var ranking = data.ranking[0];
            var patientIDs = data.patientIDs;
            var uids = data.uids;
            var numberOfClusteredColumns = data.numberOfClusteredColumns[0];
            var numberOfClusteredRows = data.numberOfClusteredRows[0];
            // var warning = data.warnings === undefined ? '' : data.warnings;
            var maxRows = data.maxRows[0];
            var geneCardsAllowed = JSON.parse(params.geneCardsAllowed);
            var rowClustering = true;
            var colClustering = true;

            var originalPatientIDs = patientIDs.slice();
            var originalUIDs = uids.slice();

            var tmpAnimationDuration = animationDuration;
            function switchAnimation(checked) {
                if (! checked) {
                    tmpAnimationDuration = animationDuration;
                    animationDuration = 0;
                } else {
                    animationDuration = tmpAnimationDuration;
                }
            }

            var gridFieldWidth = 20;
            var gridFieldHeight = 20;
            var dendrogramHeight = 300;
            var histogramHeight = 200;
            var legendWidth = 200;
            var legendHeight = 40;
            var buttonWidth = 200;
            var buttonHeight = 40;
            var buttonPadding = 20;

            var margin = {
                top: gridFieldHeight * 2 + 100 + features.length * gridFieldHeight / 2 + dendrogramHeight,
                right: gridFieldWidth + 300 + dendrogramHeight,
                bottom: 10,
                left: histogramHeight + 250
            };

            var width = gridFieldWidth * patientIDs.length;
            var height = gridFieldHeight * uids.length;

            var selectedPatientIDs = [];

            var scale = null;
            var histogramScale = null;

            function setScales() {
                scale = d3.scale.linear()
                    .domain(d3.extent(significanceValues))
                    .range((ranking === 'pval' || ranking === 'adjpval') ? [histogramHeight, 0] : [0, histogramHeight]);

                histogramScale = function (value) {
                    return (ranking === 'ttest' || ranking === 'logfold') ? scale(Math.abs(value)) : scale(value);
                };
            }
            setScales();

            function getInternalSortValue(value) {
                switch (ranking) {
                    case 'pval':
                    case 'adjpval':
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
                    .attr('width', margin.left + margin.right + (gridFieldWidth * patientIDs.length))
                    .attr('height', margin.top + margin.bottom + (gridFieldHeight * uids.length));
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
            var patientIDItems = heatmap.append('g');
            var rowSortItems = heatmap.append('g');
            var significanceSortItems = heatmap.append('g');
            var labelItems = heatmap.append('g');
            var barItems = heatmap.append('g');
            var legendItems = heatmap.append('g');
            // var warningDiv = $('#heim-heatmap-warnings').append('strong')
            //     .text(warning);

            // this code is needed for the legend generation
            var zScores = fields.map(function(d) { return d.ZSCORE; });
            var maxZScore = Math.max.apply(null, zScores);
            var minZScore = Math.min.apply(null, zScores);
            var steps = [];
            for (var i = minZScore; i < maxZScore; i+= (maxZScore - minZScore) / 200) {
                steps.push(i);
            }

            function updateHeatmap() {
                var square = squareItems.selectAll('.square')
                    .data(fields, function (d) {
                        return 'patientID-' + d.PATIENTID + '-uid-' + d.UID;
                    });

                square.enter()
                    .append('rect')
                    .attr('class', function (d) {
                        return 'square patientID-' + smartRUtils.makeSafeForCSS(d.PATIENTID) +
                            ' uid-' + smartRUtils.makeSafeForCSS(d.UID);
                    })
                    .attr('x', function (d) {
                        return patientIDs.indexOf(d.PATIENTID) * gridFieldWidth;
                    })
                    .attr('y', function (d) {
                        return uids.indexOf(d.UID) * gridFieldHeight;
                    })
                    .attr('width', gridFieldWidth)
                    .attr('height', gridFieldHeight)
                    .attr('rx', 0)
                    .attr('ry', 0)
                    .style('fill', 'white')
                    .on('mouseover', function (d) {
                        d3.select('.patientID.patientID-' + smartRUtils.makeSafeForCSS(d.PATIENTID))
                            .classed('highlight', true);
                        d3.select('.uid.uid-' + smartRUtils.makeSafeForCSS(d.UID))
                            .classed('highlight', true);

                        var html = '';
                        for (var key in d) {
                            if (d.hasOwnProperty(key)) {
                                html += key + ': ' + d[key] + '<br/>';
                            }
                        }
                        tip.show(html);
                    })
                    .on('mouseout', function () {
                        d3.selectAll('.patientID').classed('highlight', false);
                        d3.selectAll('.uid').classed('highlight', false);
                        tip.hide();
                    });


                square.transition()
                    .duration(animationDuration)
                    .attr('x', function (d) {
                        return patientIDs.indexOf(d.PATIENTID) * gridFieldWidth;
                    })
                    .attr('y', function (d) {
                        return uids.indexOf(d.UID) * gridFieldHeight;
                    })
                    .attr('width', gridFieldWidth)
                    .attr('height', gridFieldHeight);

                var colSortText = colSortItems.selectAll('.colSortText')
                    .data(patientIDs, function (d) {
                        return d;
                    });

                colSortText.enter()
                    .append('text')
                    .attr('class', 'text colSortText')
                    .attr('x', function (d, i) {
                        return i * gridFieldWidth + 0.5 * gridFieldWidth;
                    })
                    .attr('y', -2 - gridFieldHeight + 0.5 * gridFieldHeight)
                    .attr('dy', '0.35em')
                    .attr('text-anchor', 'middle')
                    .text('↑↓');

                colSortText.transition()
                    .duration(animationDuration)
                    .attr('x', function (d, i) {
                        return i * gridFieldWidth + 0.5 * gridFieldWidth;
                    })
                    .attr('y', -2 - gridFieldHeight + 0.5 * gridFieldHeight);

                var colSortBox = colSortItems.selectAll('.colSortBox')
                    .data(patientIDs, function (d) {
                        return d;
                    });

                function getValueForSquareSorting(patientID, uid) {
                    var square = d3.select('.square' + '.patientID-' + smartRUtils.makeSafeForCSS(patientID) +
                        '.uid-' + smartRUtils.makeSafeForCSS(uid));
                    return square[0][0] ? square.property('__data__').ZSCORE : Number.NEGATIVE_INFINITY;
                }

                function isSorted(arr) {
                    return arr.every(function (d, i) {
                        return i === arr.length - 1 || arr[i][1] >= arr[i + 1][1];
                    });
                }

                colSortBox.enter()
                    .append('rect')
                    .attr('class', 'box colSortBox')
                    .attr('x', function (d, i) {
                        return i * gridFieldWidth;
                    })
                    .attr('y', -2 - gridFieldHeight)
                    .attr('width', gridFieldWidth)
                    .attr('height', gridFieldHeight)
                    .on('click', function (patientID) {
                        d3.selectAll('.box').classed('sortedBy', false);
                        d3.select(this).classed('sortedBy', true);

                        var rowValues = uids.map(function (uid, idx) {
                            return [idx, getValueForSquareSorting(patientID, uid)];
                        });
                        if (isSorted(rowValues)) {
                            rowValues.sort(function (a, b) {
                                return a[1] - b[1];
                            });
                        } else {
                            rowValues.sort(function (a, b) {
                                return b[1] - a[1];
                            });
                        }
                        var sortValues = rowValues.map(function (rowValue) {
                            return rowValue[0];
                        });
                        updateRowOrder(sortValues);
                    });

                colSortBox.transition()
                    .duration(animationDuration)
                    .attr('x', function (d, i) {
                        return i * gridFieldWidth;
                    })
                    .attr('y', -2 - gridFieldHeight)
                    .attr('width', gridFieldWidth)
                    .attr('height', gridFieldHeight);

                var rowSortText = rowSortItems.selectAll('.rowSortText')
                    .data(uids, function (d) {
                        return d;
                    });

                rowSortText.enter()
                    .append('text')
                    .attr('class', 'text rowSortText')
                    .attr('transform', function (d, i) {
                        return 'translate(' + (width + 2 + 0.5 * gridFieldWidth) + ',0)' + 'translate(0,' +
                            (i * gridFieldHeight + 0.5 * gridFieldHeight) + ')rotate(-90)';
                    })
                    .attr('dy', '0.35em')
                    .attr('text-anchor', 'middle')
                    .text('↑↓');

                rowSortText.transition()
                    .duration(animationDuration)
                    .attr('transform', function (d, i) {
                        return 'translate(' + (width + 2 + 0.5 * gridFieldWidth) + ',0)' + 'translate(0,' +
                            (i * gridFieldHeight + 0.5 * gridFieldHeight) + ')rotate(-90)';
                    });

                var rowSortBox = rowSortItems.selectAll('.rowSortBox')
                    .data(uids, function (d) {
                        return d;
                    });

                rowSortBox.enter()
                    .append('rect')
                    .attr('class', 'box rowSortBox')
                    .attr('x', width + 2)
                    .attr('y', function (d, i) {
                        return i * gridFieldHeight;
                    })
                    .attr('width', gridFieldWidth)
                    .attr('height', gridFieldHeight)
                    .on('click', function (uid) {
                        d3.selectAll('.box').classed('sortedBy', false);
                        d3.select(this).classed('sortedBy', true);

                        var colValues = patientIDs.map(function (patientID, idx) {
                            return [idx, getValueForSquareSorting(patientID, uid)];
                        });
                        if (isSorted(colValues)) {
                            colValues.sort(function (a, b) {
                                return a[1] - b[1];
                            });
                        } else {
                            colValues.sort(function (a, b) {
                                return b[1] - a[1];
                            });
                        }
                        var sortValues = colValues.map(function (colValue) {
                            return colValue[0];
                        });
                        updateColOrder(sortValues);
                    });

                rowSortBox.transition()
                    .duration(animationDuration)
                    .attr('x', width + 2)
                    .attr('y', function (d, i) {
                        return i * gridFieldHeight;
                    })
                    .attr('width', gridFieldWidth)
                    .attr('height', gridFieldHeight);

                var significanceSortText = significanceSortItems.selectAll('.significanceSortText')
                    .data(['something'], function (d) {
                        return d;
                    });

                significanceSortText.enter()
                    .append('text')
                    .attr('class', 'text significanceSortText')
                    .attr('x', -gridFieldWidth - 10 + 0.5 * gridFieldWidth)
                    .attr('y', -2 - gridFieldHeight + 0.5 * gridFieldHeight)
                    .attr('dy', '0.35em')
                    .attr('text-anchor', 'middle')
                    .text('↑↓');

                significanceSortText.transition()
                    .duration(animationDuration)
                    .attr('x', -gridFieldWidth - 10 + 0.5 * gridFieldWidth)
                    .attr('y', -2 - gridFieldHeight + 0.5 * gridFieldHeight);

                var significanceSortBox = significanceSortItems.selectAll('.significanceSortBox')
                    .data(['something'], function (d) {
                        return d;
                    });

                significanceSortBox.enter()
                    .append('rect')
                    .attr('class', 'box significanceSortBox')
                    .attr('x', -gridFieldWidth - 10)
                    .attr('y', -2 - gridFieldHeight)
                    .attr('width', gridFieldWidth)
                    .attr('height', gridFieldHeight)
                    .on('click', function () {
                        d3.selectAll('.box').classed('sortedBy', false);
                        d3.select(this).classed('sortedBy', true);

                        var rowValues = significanceValues.map(function (significanceValue, idx) {
                            return [idx, getInternalSortValue(significanceValue)];
                        });

                        if (isSorted(rowValues)) {
                            rowValues.sort(function (a, b) {
                                return a[1] - b[1];
                            });
                        } else {
                            rowValues.sort(function (a, b) {
                                return b[1] - a[1];
                            });
                        }
                        var sortValues = rowValues.map(function (rowValue) {
                            return rowValue[0];
                        });
                        updateRowOrder(sortValues);
                    });

                significanceSortBox.transition()
                    .duration(animationDuration)
                    .attr('x', -gridFieldWidth - 10)
                    .attr('y', -2 - gridFieldHeight)
                    .attr('width', gridFieldWidth)
                    .attr('height', gridFieldHeight);

                var selectText = selectItems.selectAll('.selectText')
                    .data(patientIDs, function (d) {
                        return d;
                    });

                selectText.enter()
                    .append('text')
                    .attr('class', function (d) {
                        return 'text selectText patientID-' + smartRUtils.makeSafeForCSS(d);
                    })
                    .attr('x', function (d, i) {
                        return i * gridFieldWidth + 0.5 * gridFieldWidth;
                    })
                    .attr('y', -2 - gridFieldHeight * 2 + 0.5 * gridFieldHeight)
                    .attr('dy', '0.35em')
                    .attr('text-anchor', 'middle')
                    .text('□');

                selectText.transition()
                    .duration(animationDuration)
                    .attr('x', function (d, i) {
                        return i * gridFieldWidth + 0.5 * gridFieldWidth;
                    })
                    .attr('y', -2 - gridFieldHeight * 2 + 0.5 * gridFieldHeight);

                var selectBox = selectItems.selectAll('.selectBox')
                    .data(patientIDs, function (d) {
                        return d;
                    });

                selectBox.enter()
                    .append('rect')
                    .attr('class', 'box selectBox')
                    .attr('x', function (d, i) {
                        return i * gridFieldWidth;
                    })
                    .attr('y', -2 - gridFieldHeight * 2)
                    .attr('width', gridFieldWidth)
                    .attr('height', gridFieldHeight)
                    .on('click', function (patientID) {
                        selectCol(patientID);
                    });

                selectBox.transition()
                    .duration(animationDuration)
                    .attr('x', function (d, i) {
                        return i * gridFieldWidth;
                    })
                    .attr('y', -2 - gridFieldHeight * 2)
                    .attr('width', gridFieldWidth)
                    .attr('height', gridFieldHeight);

                var patientID = patientIDItems.selectAll('.patientID')
                    .data(patientIDs, function (d) {
                        return d;
                    });

                patientID.enter()
                    .append('text')
                    .attr('class', function (d) {
                        return 'patientID patientID-' + smartRUtils.makeSafeForCSS(d);
                    })
                    .attr('transform', function (d) {
                        return 'translate(' + (patientIDs.indexOf(d) * gridFieldWidth) + ',0)' +
                            'translate(' + (gridFieldWidth / 2) + ',' + (-4 - gridFieldHeight * 2) + ')rotate(-45)';
                    })
                    .style('text-anchor', 'start')
                    .text(function (d) {
                        return d;
                    });

                patientID.transition()
                    .duration(animationDuration)
                    .attr('transform', function (d) {
                        return 'translate(' + (patientIDs.indexOf(d) * gridFieldWidth) + ',0)' +
                            'translate(' + (gridFieldWidth / 2) + ',' + (-4 - gridFieldHeight * 2) + ')rotate(-45)';
                    });

                var uid = labelItems.selectAll('.uid')
                    .data(uids, function (d) {
                        return d;
                    });

                uid.enter()
                    .append('text')
                    .attr('class', function (d) {
                        return 'uid text uid-' + smartRUtils.makeSafeForCSS(d);
                    })
                    .attr('x', width + gridFieldWidth + 7)
                    .attr('y', function (d) {
                        return uids.indexOf(d) * gridFieldHeight + 0.5 * gridFieldHeight;
                    })
                    .attr('dy', '0.35em')
                    .style('text-anchor', 'start')
                    .text(function (d) {
                        return d;
                    })
                    .on('click', function(d) {
                        var genes = d.split('--');
                        genes.shift();
                        if (geneCardsAllowed) {
                            genes.forEach(function(gene) {
                                var url = 'http://www.genecards.org/cgi-bin/carddisp.pl?gene=' + gene;
                                window.open(url);
                            });
                        } else {
                            alert('To use this function you have to accept the GeneCards Terms of Use and re-create the plot');
                        }
                    });

                uid.transition()
                    .duration(animationDuration)
                    .attr('x', width + gridFieldWidth + 7)
                    .attr('y', function (d) {
                        return uids.indexOf(d) * gridFieldHeight + 0.5 * gridFieldHeight;
                    });

                var significanceIndexMap = $.map(significanceValues, function (d, i) {
                    return {significance: d, idx: i};
                });

                var bar = barItems.selectAll('.bar')
                    .data(significanceIndexMap, function (d) {
                        return d.idx;
                    });

                bar.enter()
                    .append('rect')
                    .attr('class', function (d) {
                        return 'bar idx-' + smartRUtils.makeSafeForCSS(d.idx);
                    })
                    .attr('width', function (d) {
                        return histogramScale(d.significance);
                    })
                    .attr('height', gridFieldHeight)
                    .attr('x', function (d) {
                        return -histogramScale(d.significance);
                    })
                    .attr('y', function (d, idx) {
                        return gridFieldHeight * idx;
                    })
                    .style('fill', function (d) {
                        return d.significance > 0 ? 'steelblue' : '#990000';
                    })
                    .on('mouseover', function (d) {
                        var html = 'Ranking (' + ranking + '): ' + d.significance;
                        tip.show(html);
                        d3.selectAll('.square.uid-' + smartRUtils.makeSafeForCSS(uids[d.idx]))
                            .classed('squareHighlighted', true);
                        d3.select('.uid.uid-' + smartRUtils.makeSafeForCSS(uids[d.idx]))
                            .classed('highlight', true);
                    })
                    .on('mouseout', function () {
                        tip.hide();
                        d3.selectAll('.square').classed('squareHighlighted', false);
                        d3.selectAll('.uid').classed('highlight', false);
                    });

                bar.transition()
                    .duration(animationDuration)
                    .attr('height', gridFieldHeight)
                    .attr('width', function (d) {
                        return histogramScale(d.significance);
                    })
                    .attr('x', function (d) {
                        return -histogramScale(d.significance);
                    })
                    .attr('y', function (d) {
                        return gridFieldHeight * d.idx;
                    })
                    .style('fill', function (d) {
                        return d.significance > 0 ? 'steelblue' : '#990000';
                    });

                var featurePosY = -gridFieldWidth * 2 - smartRUtils.getMaxWidth(d3.selectAll('.patientID')) -
                    features.length * gridFieldWidth / 2 - 20;

                var extraSquare = featureItems.selectAll('.extraSquare')
                    .data(extraFields, function (d) {
                        return 'patientID-' + smartRUtils.makeSafeForCSS(d.PATIENTID) +
                            '-feature-' + smartRUtils.makeSafeForCSS(d.FEATURE);
                    });

                extraSquare.enter()
                    .append('rect')
                    .attr('class', function (d) {
                        return 'extraSquare patientID-' + smartRUtils.makeSafeForCSS(d.PATIENTID) +
                            ' feature-' + smartRUtils.makeSafeForCSS(d.FEATURE);
                    })
                    .attr('x', function (d) {
                        return patientIDs.indexOf(d.PATIENTID) * gridFieldWidth;
                    })
                    .attr('y', function (d) {
                        return featurePosY + features.indexOf(d.FEATURE) * gridFieldHeight / 2;
                    })
                    .attr('width', gridFieldWidth)
                    .attr('height', gridFieldHeight / 2)
                    .attr('rx', 0)
                    .attr('ry', 0)
                    .style('fill', 'white')
                    .on('mouseover', function (d) {
                        d3.select('.patientID.patientID-' + smartRUtils.makeSafeForCSS(d.PATIENTID))
                            .classed('highlight', true);
                        d3.select('.feature.feature-' + smartRUtils.makeSafeForCSS(d.FEATURE))
                            .classed('highlight', true);
                        var html = '';
                        for (var key in d) {
                            if (d.hasOwnProperty(key)) {
                                html += key + ': ' + d[key] + '<br/>';
                            }
                        }
                        tip.show(html);
                    })
                    .on('mouseout', function () {
                        tip.hide();
                        d3.selectAll('.patientID').classed('highlight', false);
                        d3.selectAll('.feature').classed('highlight', false);
                    });

                extraSquare.transition()
                    .duration(animationDuration)
                    .attr('x', function (d) {
                        return patientIDs.indexOf(d.PATIENTID) * gridFieldWidth;
                    })
                    .attr('y', function (d) {
                        return featurePosY + features.indexOf(d.FEATURE) * gridFieldHeight / 2;
                    })
                    .attr('width', gridFieldWidth)
                    .attr('height', gridFieldHeight / 2);

                var feature = featureItems.selectAll('.feature')
                    .data(features, function (d) {
                        return d;
                    });

                feature.enter()
                    .append('text')
                    .attr('class', function (d) {
                        return 'feature text feature-' + smartRUtils.makeSafeForCSS(d);
                    })
                    .attr('x', width + gridFieldWidth + 7)
                    .attr('y', function (d) {
                        return featurePosY + features.indexOf(d) * gridFieldHeight / 2 + gridFieldHeight / 4;
                    })
                    .attr('dy', '0.35em')
                    .style('text-anchor', 'start')
                    .text(function (d) {
                        return d;
                    });

                feature.transition()
                    .duration(animationDuration)
                    .attr('x', width + gridFieldWidth + 7)
                    .attr('y', function (d) {
                        return featurePosY + features.indexOf(d) * gridFieldHeight / 2 + gridFieldHeight / 4;
                    });

                var featureSortText = featureItems.selectAll('.featureSortText')
                    .data(features, function (d) {
                        return d;
                    });

                featureSortText.enter()
                    .append('text')
                    .attr('class', 'text featureSortText')
                    .attr('transform', function (d) {
                        return 'translate(' + (width + 2 + 0.5 * gridFieldWidth) + ',0)' + 'translate(0,' +
                            (featurePosY + features.indexOf(d) * gridFieldHeight / 2 + gridFieldHeight / 4) +
                            ')rotate(-90)';
                    })
                    .attr('dy', '0.35em')
                    .attr('text-anchor', 'middle')
                    .text('↑↓');

                featureSortText.transition()
                    .duration(animationDuration)
                    .attr('transform', function (d) {
                        return 'translate(' + (width + 2 + 0.5 * gridFieldWidth) + ',0)' + 'translate(0,' +
                            (featurePosY + features.indexOf(d) * gridFieldHeight / 2 + gridFieldHeight / 4) +
                            ')rotate(-90)';
                    });

                var featureSortBox = featureItems.selectAll('.featureSortBox')
                    .data(features, function (d) {
                        return d;
                    });

                featureSortBox.enter()
                    .append('rect')
                    .attr('class', 'box featureSortBox')
                    .attr('x', width + 2)
                    .attr('y', function (d) {
                        return featurePosY + features.indexOf(d) * gridFieldHeight / 2;
                    })
                    .attr('width', gridFieldWidth)
                    .attr('height', gridFieldHeight / 2)
                    .on('click', function (feature) {
                        d3.selectAll('.box').classed('sortedBy', false);
                        d3.select(this).classed('sortedBy', true);

                        var featureValues = [];
                        var missingValues = false;
                        for (var i = 0; i < patientIDs.length; i++) {
                            var patientID = patientIDs[i];
                            var value = (-Math.pow(2, 32)).toString();
                            try {
                                var square = d3.select('.extraSquare' + '.patientID-' +
                                    smartRUtils.makeSafeForCSS(patientID) + '.feature-' +
                                    smartRUtils.makeSafeForCSS(feature));
                                value = square.property('__data__').VALUE;
                            } catch (err) {
                                missingValues = true;
                            }
                            featureValues.push([i, value]);
                        }
                        if (isSorted(featureValues)) {
                            featureValues.sort(function (a, b) {
                                var diff = a[1] - b[1];
                                return isNaN(diff) ? a[1].localeCompare(b[1]) : diff;
                            });
                        } else {
                            featureValues.sort(function (a, b) {
                                var diff = b[1] - a[1];
                                return isNaN(diff) ? b[1].localeCompare(a[1]) : diff;
                            });
                        }
                        var sortValues = [];
                        for (i = 0; i < featureValues.length; i++) {
                            sortValues.push(featureValues[i][0]);
                        }
                        if (missingValues) {
                            alert('Feature is missing for one or more patients.\n' +
                                'Every missing value will be set to lowest possible value for sorting;');
                        }
                        updateColOrder(sortValues);
                    });


                featureSortBox.transition()
                    .duration(animationDuration)
                    .attr('x', width + 2)
                    .attr('y', function (d) {
                        return featurePosY + features.indexOf(d) * gridFieldHeight / 2;
                    })
                    .attr('width', gridFieldWidth)
                    .attr('height', gridFieldHeight / 2);
            }

            function zoom(zoomLevel) {
                zoomLevel /= 100;
                d3.selectAll('.patientID')
                    .style('font-size', Math.ceil(14 * zoomLevel) + 'px');
                d3.selectAll('.selectText')
                    .style('font-size', Math.ceil(16 * zoomLevel) + 'px');
                d3.selectAll('.uid')
                    .style('font-size', Math.ceil(12 * zoomLevel) + 'px');
                d3.selectAll('.feature')
                    .style('font-size', Math.ceil(10 * zoomLevel) + 'px');
                d3.selectAll('.significanceSortText, .rowSortText, .colSortText')
                    .style('font-size', Math.ceil(14 * zoomLevel) + 'px');
                d3.selectAll('.featureSortText')
                    .style('font-size', Math.ceil(10 * zoomLevel) + 'px');
                gridFieldWidth = 20 * zoomLevel;
                gridFieldHeight = 20 * zoomLevel;
                width = gridFieldWidth * patientIDs.length;
                height = gridFieldHeight * uids.length;
                heatmap
                    .attr('width', width + margin.left + margin.right)
                    .attr('height', width + margin.top + margin.bottom);
                var temp = animationDuration;
                animationDuration = 0;
                updateHeatmap();
                reloadDendrograms();
                animationDuration = temp;
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
                d3.selectAll('.bar')
                    .map(function (d) {
                        return [d.idx, histogramScale(d.significance)];
                    }) // This line is a bit hacky
                    .sort(function (a, b) {
                        return a[1] - b[1];
                    })
                    .filter(function (d, i) {
                        return i < cutoff;
                    })
                    .each(function (d) {
                        d3.select('.bar.idx-' + smartRUtils.makeSafeForCSS(d[0]))
                            .classed('cuttoffHighlight', true);
                        d3.selectAll('.square.uid-' + smartRUtils.makeSafeForCSS(uids[d[0]]))
                            .classed('cuttoffHighlight', true);
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
                params.max_row = maxRows - cutoffLevel - 1;
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

            function selectCol(patientID) {
                var colSquares = d3.selectAll('.square.patientID-' + smartRUtils.makeSafeForCSS(patientID));
                if (colSquares.classed('selected')) {
                    var index = selectedPatientIDs.indexOf(patientID);
                    selectedPatientIDs.splice(index, 1);
                    colSquares
                        .classed('selected', false);
                    d3.select('.selectText.patientID-' + smartRUtils.makeSafeForCSS(patientID))
                        .text('□');
                } else {
                    selectedPatientIDs.push(patientID);
                    colSquares.classed('selected', true);
                    d3.select('.selectText.patientID-' + smartRUtils.makeSafeForCSS(patientID))
                        .text('■');
                }
                if (selectedPatientIDs.length !== 0) {
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
                    .duration(animationDuration)
                    .style('fill', function (d) {
                        return colorScale(1 / (1 + Math.pow(Math.E, -d.ZSCORE)));
                    });

                var featureColorSetBinary = ['#FF8000', '#FFFF00'];
                var featureColorSetSequential = [
                    'rgb(247,252,253)', 'rgb(224,236,244)', 'rgb(191,211,230)',
                    'rgb(158,188,218)', 'rgb(140,150,198)', 'rgb(140,107,177)',
                    'rgb(136,65,157)', 'rgb(129,15,124)', 'rgb(77,0,75)'
                ];
                var featureColorCategorical = d3.scale.category10();

                features.forEach(function(feature) {
                    d3.selectAll('.extraSquare.feature-' + smartRUtils.makeSafeForCSS(feature))
                        .style('fill', function (d) {
                            switch (d.TYPE) {
                                case 'binary':
                                    return featureColorSetBinary[d.VALUE];
                                case 'subset':
                                    return featureColorSetBinary[d.VALUE - 1];
                                case 'numerical':
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
                    .attr('x', function(d, i) {
                        return 2 - margin.left + buttonPadding + buttonWidth + i * legendElementWidth;
                    })
                    .attr('y', 8 - margin.top + buttonHeight * 4 + buttonPadding * 4)
                    .attr('width', Math.ceil(legendElementWidth))
                    .attr('height', legendElementHeight)
                    .style('fill', function(d) { return colorScale(1 / (1 + Math.pow(Math.E, -d))); });

                legendColor.transition()
                    .duration(animationDuration)
                    .style('fill', function(d) { return colorScale(1 / (1 + Math.pow(Math.E, -d))); });

                var legendText = legendItems.selectAll('.legendText')
                    .data(steps, function(d) { return d; });

                legendText.enter()
                    .append('text')
                    .attr('class', 'legendText')
                    .attr('x', function(d, i) {
                        return 2 - margin.left + buttonPadding + buttonWidth + i * legendElementWidth;
                    })
                    .attr('y', 8 - margin.top + buttonHeight * 4 + buttonPadding * 4 + legendHeight + 10)
                    .attr('text-anchor', 'middle')
                    .text(function(d, i) {
                        if (i === 0 || i === steps.length - 1) {
                            return Number((steps.min()).toFixed(1));
                        } else {
                            return null;
                        }
                    });

                legendText.transition()
                    .text(function(d, i) {
                        if (i === 0) {
                            return Number((steps.min()).toFixed(1));
                        } else if (i === steps.length - 1) {
                            return Number((steps.max()).toFixed(1));
                        } else {
                            return null;
                        }
                    });
            }

            function unselectAll() {
                d3.selectAll('.selectText')
                    .text('□');
                d3.selectAll('.square')
                    .classed('selected', false)
                    .attr('opacity', 1);
                selectedPatientIDs = [];
            }

            var colDendrogramVisible = false;
            var colDendrogram;

            function createColDendrogram() {
                var w = 200;
                var colDendrogramWidth = gridFieldWidth * numberOfClusteredColumns;
                var spacing = gridFieldWidth * 2 + smartRUtils.getMaxWidth(d3.selectAll('.patientID')) +
                    features.length * gridFieldHeight / 2 + 40;

                var cluster = d3.layout.cluster()
                    .size([colDendrogramWidth, w])
                    .separation(function() {
                        return 1;
                    });

                var diagonal = d3.svg.diagonal()
                    .projection(function(d) {
                        return [d.x, -spacing - w + d.y];
                    });

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
                    .attr('transform', function (d) {
                        return 'translate(' + d.x + ',' + (-spacing - w + d.y) + ')';
                    }).on('click', function (d) {
                        var previousSelection = selectedPatientIDs.slice();
                        unselectAll();
                        var leafs = d.index.split(' ');
                        for (var i = 0; i < leafs.length; i++) {
                            var patientID = patientIDs[leafs[i]];
                            selectCol(patientID);
                        }
                        if (previousSelection.sort().toString() === selectedPatientIDs.sort().toString()) {
                            unselectAll();
                        }
                    })
                    .on('mouseover', function (d) {
                        tip.show('Height: ' + d.height);
                    })
                    .on('mouseout', function () {
                        tip.hide();
                    });
                colDendrogramVisible = true;
            }

            var rowDendrogramVisible = false;
            var rowDendrogram;

            function createRowDendrogram() {
                var h = 280;
                var rowDendrogramHeight = gridFieldWidth * numberOfClusteredRows;
                var spacing = gridFieldWidth + smartRUtils.getMaxWidth(d3.selectAll('.uid')) + 20;

                var cluster = d3.layout.cluster()
                    .size([rowDendrogramHeight, h])
                    .separation(function () {
                        return 1;
                    });

                var diagonal = d3.svg.diagonal()
                    .projection(function (d) {
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
                    .attr('transform', function (d) {
                        return 'translate(' + (width + spacing + h - d.y) + ',' + d.x + ')';
                    }).on('click', function (d) {
                        var leafs = d.index.split(' ');
                        var genes = [];
                        leafs.each(function (leaf) {
                            var uid = uids[leaf];
                            var split = uid.split("--");
                            split.shift();
                            split.each(function (gene) {
                                genes.push(gene);
                            });
                        });
                        $.ajax({
                            url: 'http://biocompendium.embl.de/cgi-bin/biocompendium.cgi',
                            type: 'POST',
                            timeout: '5000',
                            async: false,
                            data: {
                                section: 'upload_gene_lists_general',
                                primary_org: 'Human',
                                background: 'whole_genome',
                                Category1: 'Human',
                                gene_list_1: 'gene_list_1',
                                SubCat1: 'hgnc_symbol',
                                attachment1: genes.join(' ')
                            }
                        }).done(function (serverAnswer) {
                            var sessionID = serverAnswer.match(/tmp_\d+/)[0];
                            var url = 'http://biocompendium.embl.de/' +
                                'cgi-bin/biocompendium.cgi?section=pathway&pos=0&background=whole_genome&session=' +
                                sessionID + '&list=gene_list_1__1&list_size=15&org=human';
                            window.open(url);
                        }).fail(function () {
                            alert('An error occurred. Maybe the external resource is unavailable.');
                        });
                    })
                    .on('mouseover', function (d) {
                        tip.show('Height: ' + d.height);
                    })
                    .on('mouseout', function () {
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

            function updateColOrder(sortValues) {
                patientIDs = sortValues.map(function (sortValue) {
                    return patientIDs[sortValue];
                });
                unselectAll();
                removeColDendrogram();
                updateHeatmap();
            }

            function updateRowOrder(sortValues) {
                var sortedUIDs = [];
                var sortedSignificanceValues = [];
                var sortedLogfoldValues = [];
                var sortedTtestValues = [];
                var sortedPvalValues = [];
                var sortedAdjpvalValues = [];
                var sortedBvalValues = [];

                sortValues.each(function (sortValue) {
                    sortedUIDs.push(uids[sortValue]);
                    sortedSignificanceValues.push(significanceValues[sortValue]);
                    sortedLogfoldValues.push(logfoldValues[sortValue]);
                    sortedTtestValues.push(ttestValues[sortValue]);
                    sortedPvalValues.push(pvalValues[sortValue]);
                    sortedAdjpvalValues.push(adjpvalValues[sortValue]);
                    sortedBvalValues.push(bvalValues[sortValue]);
                });
                uids = sortedUIDs;
                significanceValues = sortedSignificanceValues;
                logfoldValues = sortedLogfoldValues;
                ttestValues = sortedTtestValues;
                pvalValues = sortedPvalValues;
                adjpvalValues = sortedAdjpvalValues;
                bvalValues = sortedBvalValues;

                removeRowDendrogram();
                updateHeatmap();
                animateCutoff();
            }

            function transformClusterOrderWRTInitialOrder(clusterOrder, initialOrder) {
                return clusterOrder.map(function (d) {
                    return initialOrder.indexOf(d);
                });
            }

            function getInitialRowOrder() {
                return uids.map(function (uid) {
                    return originalUIDs.indexOf(uid);
                });
            }

            function getInitialColOrder() {
                return patientIDs.map(function (patientID) {
                    return originalPatientIDs.indexOf(patientID);
                });
            }

            var lastUsedClustering = null;

            function cluster(clustering) {
                if (!lastUsedClustering && typeof clustering === 'undefined') {
                    return; // Nothing should be done if clustering switches are turned on without clustering type set.
                }
                d3.selectAll('.box').classed('sortedBy', false);
                clustering = (typeof clustering === 'undefined') ? lastUsedClustering : clustering;
                var clusterData = data[clustering];
                if (rowClustering && numberOfClusteredRows > 0) {
                    rowDendrogram = JSON.parse(clusterData[3]);
                    updateRowOrder(transformClusterOrderWRTInitialOrder(clusterData[1], getInitialRowOrder()));
                    createRowDendrogram(rowDendrogram);
                } else {
                    removeRowDendrogram();
                }
                if (colClustering && numberOfClusteredColumns > 0) {
                    colDendrogram = JSON.parse(clusterData[2]);
                    updateColOrder(transformClusterOrderWRTInitialOrder(clusterData[0], getInitialColOrder()));
                    createColDendrogram(colDendrogram);
                } else {
                    removeColDendrogram();
                }
                lastUsedClustering = clustering;
            }

            function switchRowClustering() {
                rowClustering = !rowClustering;
                cluster();
            }

            function switchColClustering() {
                colClustering = !colClustering;
                cluster();
            }

            function changeRanking(method) {
                ranking = method;
                switch(method) {
                    case 'pval':
                        significanceValues = pvalValues;
                        break;
                    case 'adjpval':
                        significanceValues = adjpvalValues;
                        break;
                    case 'ttest':
                        significanceValues = ttestValues;
                        break;
                    case 'logfold':
                        significanceValues = logfoldValues;
                        break;
                    default:
                        significanceValues = bvalValues;
                }
                setScales();
                updateHeatmap();
            }

            function init() {
                updateHeatmap();
                reloadDendrograms();
                updateColors('redGreen');
            }

            init();

            controlElements.createD3Switch({
                location: heatmap,
                onlabel: 'Animation ON',
                offlabel: 'Animation OFF',
                x: 2 - margin.left,
                y: 8 - margin.top,
                width: buttonWidth,
                height: buttonHeight,
                callback: switchAnimation,
                checked: true
            });

            controlElements.createD3Slider({
                location: heatmap,
                label: 'Zoom in %',
                x: 2 - margin.left + buttonPadding + buttonWidth,
                y: 8 - margin.top - 10,
                width: buttonWidth,
                height: buttonHeight,
                min: 1,
                max: 200,
                init: 100,
                callback: zoom,
                trigger: 'dragend'
            });

            controlElements.createD3Button({
                location: heatmap,
                label: 'Apply Cutoff',
                x: 2 - margin.left,
                y: 8 - margin.top + buttonHeight + buttonPadding,
                width: buttonWidth,
                height: buttonHeight,
                callback: cutoff
            });

            controlElements.createD3Slider({
                location: heatmap,
                label: 'Cutoff',
                x: 2 - margin.left + buttonPadding + buttonWidth,
                y: 8 - margin.top + buttonHeight + buttonPadding - 10,
                width: buttonWidth,
                height: buttonHeight,
                min: 0,
                max: maxRows - 2,
                init: 0,
                callback: animateCutoff,
                trigger: 'dragend'
            });

            controlElements.createD3Switch({
                location: heatmap,
                onlabel: 'Clustering rows ON',
                offlabel: 'Clustering rows OFF',
                x: 2 - margin.left,
                y: 8 - margin.top + buttonHeight * 3 + buttonPadding * 3,
                width: buttonWidth,
                height: buttonHeight,
                callback: switchRowClustering,
                checked: rowClustering
            });

            controlElements.createD3Switch({
                location: heatmap,
                onlabel: 'Clustering columns ON',
                offlabel: 'Clustering columns OFF',
                x: 2 - margin.left + buttonPadding + buttonWidth,
                y: 8 - margin.top + buttonHeight * 3 + buttonPadding * 3,
                width: buttonWidth,
                height: buttonHeight,
                callback: switchColClustering,
                checked: colClustering
            });

            if (ranking === 'bval' ||
                ranking === 'pval' ||
                ranking === 'adjpval' ||
                ranking === 'logfold' ||
                ranking === 'ttest') {
                controlElements.createD3Dropdown({
                    location: heatmap,
                    label: 'Ranking Method',
                    x: 2 - margin.left,
                    y: 8 - margin.top + buttonHeight * 4 + buttonPadding * 4,
                    width: buttonWidth,
                    height: buttonHeight,
                    items: [
                        {
                            callback: function () {
                                changeRanking('bval');
                            },
                            label: 'B Value'
                        },
                        {
                            callback: function () {
                                changeRanking('ttest');
                            },
                            label: 'T Test'
                        },
                        {
                            callback: function () {
                                changeRanking('logfold');
                            },
                            label: 'Logfold'
                        },
                        {
                            callback: function () {
                                changeRanking('pval');
                            },
                            label: 'p Value'
                        },
                        {
                            callback: function () {
                                changeRanking('adjpval');
                            },
                            label: 'adj. p Value'
                        }
                    ]
                });
            }


            controlElements.createD3Dropdown({
                location: heatmap,
                label: 'Heatmap Coloring',
                x: 2 - margin.left,
                y: 8 - margin.top + buttonHeight * 2 + buttonPadding * 2,
                width: buttonWidth,
                height: buttonHeight,
                items: [
                    {
                        callback: function () {
                            updateColors('redGreen');
                        },
                        label: 'Red to Green Schema'
                    },
                    {
                        callback: function () {
                            updateColors('redBlue');
                        },
                        label: 'Red to Blue Schema'
                    },
                    {
                        callback: function () {
                            updateColors('blueScale');
                        },
                        label: 'Blue Schema'
                    },
                    {
                        callback: function () {
                            updateColors('greenScale');
                        },
                        label: 'Green Schema'
                    }
                ]
            });

            controlElements.createD3Dropdown({
                location: heatmap,
                label: 'Heatmap Clustering',
                x: 2 - margin.left + buttonPadding + buttonWidth,
                y: 8 - margin.top + buttonHeight * 2 + buttonPadding * 2,
                width: buttonWidth,
                height: buttonHeight,
                items: [
                    {
                        callback: function () {
                            cluster('hclustEuclideanAverage');
                        },
                        label: 'Hierarch.-Eucl.-Average'
                    },
                    {
                        callback: function () {
                            cluster('hclustEuclideanComplete');
                        },
                        label: 'Hierarch.-Eucl.-Complete'
                    },
                    {
                        callback: function () {
                            cluster('hclustEuclideanSingle');
                        },
                        label: 'Hierarch.-Eucl.-Single'
                    },
                    {
                        callback: function () {
                            cluster('hclustManhattanAverage');
                        },
                        label: 'Hierarch.-Manhat.-Average'
                    },
                    {
                        callback: function () {
                            cluster('hclustManhattanComplete');
                        },
                        label: 'Hierarch.-Manhat.-Complete'
                    },
                    {
                        callback: function () {
                            cluster('hclustManhattanSingle');
                        },
                        label: 'Hierarch.-Manhat.-Single'
                    }
                ]
            });
        }

    }]);
