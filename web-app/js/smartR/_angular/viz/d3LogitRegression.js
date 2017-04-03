//# sourceURL=d3LogitRegression.js

"use strict";

window.smartRApp.directive("logitRegressionPlot", [
    "smartRUtils",
    "rServeService",
    function(smartRUtils, rServeService) {

        return {
            restrict: "E",
            scope: {
                data: "=",
                width: "@",
                height: "@"
            },
            link: function (scope, element) {
                scope.$watch("data", function() {
                    $(element[0]).empty();
                    if (! $.isEmptyObject(scope.data)) {
                        smartRUtils.prepareWindowSize(scope.width, scope.height);
                        createLogitViz(scope, element[0]);
                    }
                });
            }
        };

	function createLogitViz(scope, root){
		// layout parameter setup
		var ticks = 10;
		var plotSpace = 32;
		var xScale = 0.8;
		var yScale = 0.65;

		var w = parseInt(scope.width, 10);
		var h = parseInt(scope.height, 10);

		var margin = {top: 64, right: 256, bottom: 128, left: 128};
		var width = w * xScale - margin.left - margin.right;
		var height = h * yScale - margin.top - margin.bottom;

		var origin = {x: plotSpace + margin.left, y: plotSpace + margin.top};
		var dimensions = {width: width, height: height, margin: margin, origin: origin};

		var raw,
			data,
			nodes,
			fitted,
			residuals,
			deviance,
			pvalue;

		var xArrLabel = scope.data.xArrLabel[0];
		var yArrLabel = scope.data.yArrLabel[0];

		data = scope.data;
		updateData(data);

		var svg = d3.select(root)
			.append("svg")
			.attr("width", width + margin.left + margin.right)
			.attr("height", height + margin.top + margin.bottom);

		var domX = d3.extent(nodes, function(d){ return d.x; });
		var domY = d3.extent(nodes, function(d){ return d.y; });
		var domains = {domX: domX, domY: domY};
		
		var axdomX = d3.extent(raw, function(d){ return d.x; });
		var axdomY = d3.extent(raw, function(d){ return d.y; });
		var axisdomains = {domX: axdomX, domY: axdomY};

		updatePlotData();
		setupVisuals();
		scatterPlot(nodes);
		updatePlots();

		function removePlot() {
			d3.select(root).selectAll("*").remove();
			d3.selectAll(".d3-tip").remove();
		}

		function cleanPlotted(){
			d3.select(root).selectAll(".sr-logit-curve").remove();
			d3.select(root).selectAll(".sr-reference-line").remove();
			d3.select(root).selectAll(".sr-distr-hist").remove();
			d3.select(root).selectAll(".legend").remove();
		}

		function updateData(data) {
			raw = data.raw;
			nodes = data.data;
			fitted = data.fitted;
			residuals = data.residuals;
			deviance = data.deviance;
			pvalue = data.pvalue[1]["Pr(>Chi)"];
		}

		function updatePlotData(){
			// add rescaled coordinates for plotting
			fitted = createPlotData(fitted);
			nodes = createPlotData(nodes);

			// data augmentation 
			for(var i=0; i<nodes.length; ++i){
				nodes[i].patientID = data.data[i].patientID;
				nodes[i].fittedY = fitted[i].ploty;
				nodes[i].residual = residuals[i]; 
				nodes[i].label = "<div align=center>" + 
					smartRUtils.shortenConcept(xArrLabel) + ": " + raw[i].x.toFixed(3) + "<br>" +
					smartRUtils.shortenConcept(yArrLabel) + ": " + raw[i].y.toFixed(3) + "<br>" +
					"Residual: " + residuals[i] + "<br>" +
					"Patient ID: " + data.data[i].patientID + "</div>";
			}
		}

		function setupVisuals(){
			setupGrid();
			updateAxis();
			setupBrush(data);
			updateLegend();
		}

		function updatePlots(){
			plotRegressionCurve(fitted);
			updateReferenceLines();
			updateDistributionHistogram(nodes);
		}

		function updateReferenceLines(){
			svg.append("g")
				.selectAll("line")
				.data(nodes)
				.enter()
				.append("line")
					.attr("class", "sr-reference-line")
					.attr("id", function(d){ 
						var id = "refline-" + Math.round(d.plotx) + "-" + Math.round(d.ploty);
						return id; 
					})
					.attr("x1", function(d){ return d.plotx; })
					.attr("x2", function(d){ return d.plotx; })
					.attr("y1", function(d){ return d.ploty; })
					.attr("y2", function(d){ return d.fittedY; })
					.style("display", "none");
		}

		function setupBrush(data){
			var margin = dimensions.margin;
			var width = dimensions.width;
			var height = dimensions.height;
			var origin = dimensions.origin;

			var brushFunction = function(){
				var extent = brush.extent();
				var x0 = extent[0][0];
				var x1 = extent[1][0];
				var y0 = extent[0][1];
				var y1 = extent[1][1];

				var all = d3.selectAll(".sr-point");

				var selected = all.filter(function(d){
					return (d.plotx >= x0) && (d.plotx <= x1) && (d.ploty >= y0) && (d.ploty <= y1);
				});

				var unselected = all.filter(function(d){
					return (d.plotx < x0) || (d.plotx > x1) || (d.ploty < y0) || (d.ploty > y1);
				});

				selected.classed("selected", true);
				unselected.classed("selected", false);

				var patientIDs = selected.data().map(function(d){ return d.patientID; });
				updateRegression(patientIDs, data, false);
			};

			var brush = d3.svg.brush()
				.x(d3.scale.identity().domain([margin.left, width + origin.x + plotSpace]))
				.y(d3.scale.identity().domain([margin.top - plotSpace, height + origin.y]))
				.on("brushend", brushFunction);

			svg.append("g")
				.attr("class", "brush")
				.on("mousedown", function() {
					return d3.event.button === 2 ? d3.event.stopImmediatePropagation() : null;
				})
				.call(brush);
		}

		// plot regression curve
		// data points are sorted before plotting
		function plotRegressionCurve(data){
			// input data not necessarily sorted
			var dataSorted = data.sort(function(a, b){ return (a.x - b.x); });
			var lineGen = d3.svg.line()
				.x(function(d){ return d.plotx; })
				.y(function(d){ return d.ploty; })
				.interpolate("cardinal");

			svg.append("g")
				.append("path")
				.attr("class", "sr-logit-curve")
				.attr("d", lineGen(dataSorted));
		}

		function updateLegend() {
			var origin = dimensions.origin;

			var residualmin = d3.min(nodes.map(function(d){ return d.residual; }));
			var residualmax = d3.max(nodes.map(function(d){ return d.residual; }));
			var residualrange = "[" + residualmin.toFixed(3) +  ", " + residualmax.toFixed(3) + "]";

			var _LEGEND_LINE_SPACING = 15;
			var _LEGEND_RECT_SIZE = 10;
			var labels = [
				"X Transformation: " + data.transformationx[0],
				"Y Transformation: " + data.transformationy[0],
				"Deviance: " + deviance,
				"p-Value: " + pvalue,
				"Residual Range: " + residualrange
			];

			if(data.forceNormalize[0] && (data.transformationy[0] === "raw")){
				labels.push("Values Normalized to Range [0, 1]");
			}

			var legend = svg.selectAll(".legend")
				.data(labels, function(d) { return d; });

			var legendEnter = legend.enter()
				.append("g")
				.attr("class", "legend")
				.attr("transform", function(d, i) {
					return "translate(" + (dimensions.width+origin.x + 40) + "," + (i * (_LEGEND_RECT_SIZE + _LEGEND_LINE_SPACING)) + ")";
				});

			legendEnter.append("text");
			legend.select("text")
				.attr("y", 9)
				.text(function(d) { return d; });

			legend.exit()
				.remove();
		}

		function updateRegression(patientIDs, data, init){
			if (! init) {
				patientIDs = patientIDs.length !== 0 ? patientIDs : d3.selectAll("sr-point").data().map(function(d) {
					return d.patientID;
				});
			}

			var args = { transformationx: data.transformationx[0], transformationy: data.transformationy[0], selectedPatientIDs: patientIDs };
			rServeService.startScriptExecution({
				taskType: "run",
				arguments: args
			}).then(
				function (response) {
					var results = JSON.parse(response.result.artifacts.value);
					if (init) {
						removePlot();
						updateData(results);
					} else {
						cleanPlotted();
						updateData(results);
						updatePlotData();

						updateLegend();
						updatePlots();
					}
				},
				function (response) {
					console.error(response);
				}
			);
		}

		// add designated points
		function scatterPlot(data){
			var tip = d3.tip()
				.attr("class", "d3-tip")
				.html(function(d){ return d; });

			svg.call(tip);
			svg.append("g")
				.selectAll("circle")
				.data(data)
				.enter()
				.append("circle")
					.attr("class", "sr-point")
					.attr("cx", function(d){ return d.plotx; })
					.attr("cy", function(d){ return d.ploty; })
					.attr("r", "5")
					.on("mouseover", function(d){
						d3.select(this)
							.classed("hover", true)
							.attr("r", "7");
						
						var id = "#refline-" + Math.round(d.plotx) + "-" + Math.round(d.ploty);
						d3.select(id)
							.style("display", "block");

						tip.direction( d.residual < 0? "s" : "n" ) 
							.offset( d.residual < 0? [3, 0] : [-3, 0] )
							.show( d.label );
					})
					.on("mouseout", function(d){
						d3.select(this)
							.classed("hover", false)
							.attr("r", "5");

						var id = "#refline-" + Math.round(d.plotx) + "-" + Math.round(d.ploty);
						d3.select(id).style("display", "none");

						tip.hide();
					});
		}

		// setup sidebar and bottombar for distribution
		function updateDistributionHistogram(data){
			var histSize = 4;
			
			var margin = dimensions.margin;
			var origin = dimensions.origin;
			var height = dimensions.height;

			var countx = [];
			var county = [];

			var filterx = function(d){ return (d.plotx === dix); };
			var filtery = function(d){ return (d.ploty === diy); };

			for(var i=0; i<data.length; ++i){
				var dix = data[i].plotx;
				var diy = data[i].ploty;
				countx.push({value: dix, count: data.filter(filterx).length });
				county.push({value: diy, count: data.filter(filtery).length });
			}
			var maxcountX = d3.max(county.map(function(d){ return d.count; }));
			var maxcountY = d3.max(county.map(function(d){ return d.count; }));

			svg.append("g")
				.selectAll("rect")
				.data(countx)
				.enter()
				.append("rect")
					.attr("class", "sr-distr-hist")
					.attr("x", function(d){ return (d.value - histSize/2); })
					.attr("y", function(){ return height+origin.y; })
					.attr("width", histSize)
					.attr("height", function(d){ return (d.count/maxcountX)*(margin.top); });

			svg.append("g")
				.selectAll("rect")
				.data(county)
				.enter()
				.append("rect")
					.attr("class", "sr-distr-hist")
					.attr("x", function(d){ return margin.left-(d.count/maxcountY)*margin.left; })
					.attr("y", function(d){ return (d.value - histSize/2); })
					.attr("width", function(d){ return (d.count/maxcountY)*margin.left; })
					.attr("height", histSize);
		}

		function setupGrid(){
			var height = dimensions.height;
			var width = dimensions.width;
			var origin = dimensions.origin;
			var margin = dimensions.margin;

			var spacingx = width / ticks;
			var spacingy = height / ticks;

			var grid = [];
			for(var i=0; i<=ticks; ++i){
				grid.push(i);
			}

			// x axis grid
			svg.append("g")
				.selectAll("line")
				.data(grid)
				.enter()
				.append("line")
				.attr("class", "sr-grid-line")
				.attr("x1", origin.x)
				.attr("x2", origin.x + width)
				.attr("y1", function(d, i){ return Math.round((i*spacingy) + margin.top); })
				.attr("y2", function(d, i){ return Math.round((i*spacingy) + margin.top); });

			// y axis
			svg.append("g")
				.selectAll("line")
				.data(grid)
				.enter()
				.append("line")
				.attr("class", "sr-grid-line")
				.attr("x1", function(d, i){ return Math.round((i*spacingx) + origin.x); })
				.attr("x2", function(d, i){ return Math.round((i*spacingx) + origin.x); })
				.attr("y1", margin.top)
				.attr("y2", height + margin.top);
		}

		// setup and add axes
		function updateAxis(){
			var width = dimensions.width;
			var height = dimensions.height;

			var margin = dimensions.margin;
			var origin = dimensions.origin;

			var xScale = d3.scale.linear()
					.domain(axisdomains.domY)
					.range([0, width]);

			var yScale;
			yScale = d3.scale.linear()
				.domain(axisdomains.domX)
				.range([0, height]);

			var xAxis = d3.svg.axis()
				.scale(xScale)
				.ticks(ticks)
				.orient("top");

			var yAxis = d3.svg.axis()
				.scale(yScale)
				.ticks(ticks)
				.orient("right");

			svg.append("g")
				.attr("class", "axis-x")
				.attr("transform", "translate("+origin.x+"," + (height + origin.y) + ")")
				.call(xAxis);

			svg.append("g")
				.attr("class", "axis-y")
				.attr("transform", "translate("+margin.left+", "+margin.top+")")
				.call(yAxis);

			var labelx = svg.append("text")
                .attr("class", "label-x")
                .text(smartRUtils.shortenConcept(xArrLabel));
			var boundsx = labelx.node().getBoundingClientRect();
			labelx.attr("transform", "translate("  + ((width / 2) - (boundsx.width/2) + origin.x) + ", " + (margin.top - plotSpace) + ")");

			var labely = svg.append("text")
                .attr("class", "label-y")
                .text(smartRUtils.shortenConcept(yArrLabel));
			var boundsy = labely.node().getBoundingClientRect();
			labely.attr("transform", "translate("  + (width + origin.x + plotSpace) + "," + ((height / 2) + margin.top - (boundsy.width/2)) + ") rotate(-90)");
		}

		// transform data to plot coordinates
		function createPlotData(data){
			var width = dimensions.width;
			var height = dimensions.height;

			var origin = dimensions.origin;
			var margin = dimensions.margin;

			var minX, maxX, minY, maxY;
			if(axisdomains){
				minX = domains.domX[0];
				maxX = domains.domX[1];

				minY = domains.domY[0];
				maxY = domains.domY[1];
			}else{
				// fallback
				minX = d3.min(data.map(function(d){ return d.x; }));
				maxX = d3.max(data.map(function(d){ return d.x; }));

				minY = d3.min(data.map(function(d){ return d.y; }));
				maxY = d3.max(data.map(function(d){ return d.y; }));
			}

			var rangeX = (maxX-minX)? maxX-minX : 1;
			var rangeY = (maxY-minY)? maxY-minY : 1;

			return data.map(function(d){
				return {
					x: d.x,
					y: d.y,
					plotx: width * ((d.x - minX) / rangeX) + origin.x,
					ploty: height - height * ((d.y - minY) / rangeY) + margin.top
				};
			});
		}
	}

}]);
