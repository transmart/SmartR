//# sourceURL=ppmiDemo.js

'use strict';

window.smartRApp.directive('ppmiDemo', [
    'smartRUtils',
    'rServeService',
    '$rootScope',
    function(smartRUtils, rServeService, $rootScope) {

        return {
            restrict: 'E',
            scope: {
                data: '=',
                show: '='
            },
            templateUrl: $rootScope.smartRPath +  '/js/smartR/_angular/templates/ppmidemo.html',
            link: function (scope, element) {
                var vizDiv1 = element.find('#sr-ppmi-s1')[0];
                var vizDiv2 = element.find('#sr-ppmi-s2')[0];

                scope.$watch('show', function () {
                    $(vizDiv1).empty();
                    $(vizDiv2).empty();
                    if (scope.show) {
                        createppmi(scope, vizDiv1, 1);
                        createppmi(scope, vizDiv2, 2);
                    }
                }, true);
            }
        };

        function linspace(a,b,n) {
            if(typeof n === "undefined") n = Math.max(Math.round(b-a)+1,1);
            if(n<2) { return n===1?[a]:[]; }
            var i,ret = Array(n);
            n--;
            for(i=n;i>=0;i--) { ret[i] = (i*b+(n-i)*a)/n; }
            return ret;
        }

        function createppmi(scope, vizDiv, subset) {
            var cf = crossfilter(scope.data);
            var bySubset = cf.dimension(function(d) { return d.subset; });
            var byChr = cf.dimension(function(d) { return d.chr; });

            bySubset.filterExact(subset);

            var chrs = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'x', 'y'];
            var boxColor = [];
            var allColors = linspace(0, 360, chrs.length);
            var data = [];
            var yValues = [];

            //Colors

            chrs.forEach(function(d, i) {
                var color = 'hsl('+ allColors[i] +',50%'+',50%)';
                boxColor.push(color);
            });

            //Create Y Values

            chrs.forEach(function(chr) {
                byChr.filterExact(chr);
                var frqs = byChr.top(Infinity).map(function(d) { return d.frq; });
                yValues.push(frqs);
                byChr.filterAll();
            });

            //Create Traces

            chrs.forEach(function(chr, i) {
                var result = {
                    y: yValues[i],
                    type:'box',
                    name: 'Chr' + chr,
                    marker:{
                        color: boxColor[i]
                    }
                };
                data.push(result);
            });

            //Format the layout

            var layout = {
                title: 'Subset ' + subset,
                xaxis: {
                    showgrid: false,
                    zeroline: false,
                    tickangle: 60,
                },
                yaxis: {
                    title: 'Allele Frq.',
                    zeroline: false,
                    gridcolor: 'white',
                    range: [0,1],
                }
            };


            Plotly.newPlot(vizDiv, data, layout);
            bySubset.filterAll();
        }
    }]);

