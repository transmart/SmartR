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

            var width = 1000;
            var height = 1000;
            var margin = {
                top: 200,
                bottom: 20,
                right: 20,
                left: 200
            };
            var boxSize = width / smartRUtils.unique(getValuesForDimension(bySubject)).length;

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
                    .attr('transform', 'translate(0,' + (height + 5) + ')rotate(45)')
                    .attr('text-anchor', 'start')
                    .text(function(d) { return d; });

                // UPDATE g
                verticalGridLine.attr('transform', function(d) {
                    return 'translate(' + ((subjects.indexOf(d) + 0.5) * boxSize) + ', 0)';
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
                    .attr('transform', 'translate(' + (- margin.left) + ', 0)')
                    .attr('text-anchor', 'end')
                    .text(function(d) { return d; });

                // UPDATE g
                horizontalGridLine.attr('transform', function(d) {
                    return 'translate(0, ' + ((genes.indexOf(d) + 0.5) * boxSize) + ')';
                });
            })();


            (function drawBoxes() {

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
