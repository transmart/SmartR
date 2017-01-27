//# sourceURL=d3Survival.js

'use strict';

window.smartRApp.directive('survivalPlot', [
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
						createSurvivalViz(scope, element[0]);
					}
				});
			}
		};
		
		function createSurvivalViz(scope, root) {
			
			/* User Settings */
			var plotWidth = parseInt(scope.data.plot_width);
			var plotHeight = parseInt(scope.data.plot_height);
			var timeOut = scope.data.time_out;
			var legendType = scope.data.legend_type;
			var legendPosition = scope.data.legend_position;
			var showRiskTable = String(scope.data.show_risk_table);
			showRiskTable = (showRiskTable == "true") ? true : false;
			var mergeSubsets = String(scope.data.merge_subsets);
			mergeSubsets = (mergeSubsets == "true") ? true : false;
			var mergeCategories = String(scope.data.merge_categories);
			mergeCategories = (mergeCategories == "true") ? true : false;
			var showLegend = true;
			
			/* Visualization Settings */
			var xLabel = smartRUtils.shortenConcept(new String(scope.data.x_label));
			var selectedSubsets = scope.data.selected_subsets;
			var selectedCategories = scope.data.selected_categories;
			
			/* Drawing Elements */
			var kaplanContainer;
			var kaplanPlot;
			var kaplanLegend;
			var kaplanRiskTableContainer;
			var kapalnRiskTableContainerBody;
			var kaplanRiskTable;
			
			function makeLegendLabels() {
				var legendLabels = [];
				if(mergeSubsets && mergeCategories) {
					showLegend = false;
					legendLabels.push(""); //BugFix for RiskTableHeight
				} else {
					if(selectedCategories.length > 0) {
						if(!mergeSubsets && !mergeCategories) {
							for(var i = 0; i < selectedSubsets.length; i++) {
								for(var j = 0; j < selectedCategories.length; j++) {
									legendLabels.push(selectedSubsets[i] + " - " + smartRUtils.shortenConcept(new String(selectedCategories[j])));
								}
							}
						}
						if(!mergeSubsets && mergeCategories) {
							for(var i = 0; i < selectedSubsets.length; i++) {
								legendLabels.push(selectedSubsets[i]);
							}
						}
						if(mergeSubsets && !mergeCategories) {
							for(var i = 0; i < selectedCategories.length; i++) {
								legendLabels.push(smartRUtils.shortenConcept(new String(selectedCategories[i])));
							}
						}
					} else {
						if(!mergeSubsets) {
							for(var i = 0; i < selectedSubsets.length; i++) {
								legendLabels.push(selectedSubsets[i]);
							}
						}
					}
				}
				
				return legendLabels;
			}
			
			/* Padding & Margin Settings */
			var kaplanContainerPadding = {
				top: 10,
				right: 10,
				bottom: 10,
				left: 10
			};
			var kaplanPlotPadding = {
				top: 20,
				right: 50,
				bottom: 20,
				left: 50
			};
			var kaplanLegendPadding = {
				top: 10,
				right: 15,
				bottom: 10,
				left: 15
			};
			var kaplanLegendMargin = {
				top: 10,
				right: 10,
				bottom: 10,
				left: 10
			};
			var kaplanLegendBoxMargin = {
				top: 5,
				right: 5,
				bottom: 5,
				left: 5
			};
			var kaplanRiskTablePadding = {
				top: 10,
				right: 0,
				bottom: 10,
				left: 0
			};
			
			var kaplanPlotLabelHeight = 15;
			
			// Width Variables
			var kaplanPlotInnerWidth = 0;
			var kaplanPlotWidth = 0;
			var kaplanLegendWidth = 0;
			var kaplanContainerWidth = 0;
			var kaplanRiskTableWidth = 0;
			var kaplanRiskTableContainerWidth = 0;
			
			// Legend Settings
			var legendLabels = makeLegendLabels();
			var kaplanLegendCount = legendLabels.length;
			var kaplanLegendBoxHeight = 10;
			var kaplanLegendBoxWidth = 10;
			var kaplanLegendLineHeight = kaplanLegendBoxHeight + kaplanLegendBoxMargin.top + kaplanLegendBoxMargin.bottom;
			
			// Height Variables
			var kaplanLegendHeight = 0;
			var kaplanPlotInnerHeight = 0;
			var kaplanPlotHeight = 0;
			var kaplanContainerHeight = 0;
			var kaplanRiskTableContainerHeight = 0;
			var kaplanRiskTableHeight = 0;
			
			// Position Variable
			// Will be set in function calculatePositions()
			var positions = {
				kaplanContainerX: 0,
				kaplanContainerY: 0,
				kaplanPlotX: 0,
				kaplanPlotY: 0,
				kaplanPlotXAxisX: 0,
				kaplanPlotXAxisY: 0,
				kaplanPlotXLabelX: 0,
				kaplanPlotXLabelY: 0,
				kaplanPlotYAxisX: 0,
				kaplanPlotYAxisY: 0,
				kaplanPlotYLabelX: 0,
				kaplanPlotYLabelY: 0,
				kaplanLegendX: 0,
				kaplanLegendY: 0,
				kaplanRiskTableX: 0,
				kaplanRiskTableY: 0,
			}
			
			/* Design */
			var colors = [
				'#e41a1c',
				'#377eb8',
				'#984ea3',
				'#ff7f00',
				'#ffff33',
				'#a65628',
				'#f781bf',
				'#4daf4a'
			];
			
			/* Data */
			var max = 0;
			var min= 10000000000000000000;
			var survivalData = scope.data.survival_data;
			
			/* Computed Data progression, survival, prob, censored  */
			for(var a = 0; a < survivalData.length; a++){
				for (var b = 0; b < survivalData[a].length; b++){
					var reed = survivalData[a][b];
					var brad = (b>0) ? survivalData[a][b-1].n - reed.d : reed.n;
					reed.progression = reed.d/reed.n;
					reed.survival = 1 - reed.progression;
					reed.prob = (b == 0) ? reed.survival : survivalData[a][b-1].prob * reed.survival;
					max = (max < reed.t) ? reed.t : max;
					min = (min < reed.t) ? min : reed.t;
					reed.censored = (reed.n < brad) ? true : false;
				}
			}
			
			/* Begin d3.js */
				
				// Define domains for the axes
				var xDomain = [min, max];
				var yDomain = [1, 0];
				
				//Scalar functions
				var x = d3.scale.linear().range([0, plotWidth]).domain(xDomain).nice();
				var y = d3.scale.linear().range([0, plotHeight]).domain(yDomain);
				
				// This chart will display years as integers, and populations with thousands separators
				var formatY = d3.format(",");
				var formatX = d3.format(".");
				
				// TickValues
				var tickCountX = x.ticks().length;
				var tickCountY = y.ticks().length;
				var maxXValue = x.domain()[1];
				
				//Define axes
				var xAxis = d3.svg.axis().scale(x).innerTickSize(-plotHeight).outerTickSize(2).tickPadding(6).orient("bottom").tickFormat(formatX);
				var yAxis = d3.svg.axis().scale(y).innerTickSize(-plotWidth).outerTickSize(2).tickPadding(6).orient('left');
				
				//This is the accessor function
				var lineFunction = d3.svg.line()
										.x(function(d) { return x(d.t); })
										.y(function(d) { return y(d.prob); })
										.interpolate("step-before");
				
				/* Drawing starts here */
					
					// Draw the kaplan container (contains the kaplan plot, legend and )
					kaplanContainer = d3.select(root).append('svg').attr('id', "kaplanContainer");
					
					// Draw the kaplan plot (contains the paths, axes and labels)
					kaplanPlot = kaplanContainer.append('svg').attr('id', "kaplanPlot");
					
					// Draw the x-axis
					var theXAxis = kaplanPlot.append("g").attr("class", "axis").call(xAxis);
						
					// Add the text label for the x axis
					var theXLabel = kaplanPlot.append("text").attr('class', 'axis-label').text(xLabel + " [" + timeOut + "]");
					
					//Draw the y-axis
					var theYAxis = kaplanPlot.append("g").attr("class", "axis").call(yAxis);
						
					// Add the text label for the Y axis
					var theYLabel = kaplanPlot.append("text").attr('class', 'axis-label').text("Fraction of Patients").attr("transform", "rotate(-90)");
					
					// Draw the tooltip-container
					var hoverDiv = d3.select(root).append("div").attr("class", "tooltip");
					
					// Draw the lines
					var lines = [];
					for(var a=0; a < survivalData.length; a++){
						var line = kaplanPlot.append("path")
												.attr("d", lineFunction(survivalData[a]))
												.attr("class", "active")
												.attr("id", "line-" + a)
												.attr("stroke", colors[a])
												.attr("stroke-width", 3)
												.attr("fill", "none")
												.attr('transform', 'translate(' + kaplanPlotPadding.left + ',' + kaplanPlotPadding.top + ')')
												.on('mouseover', function () {
													//on mouseover of each line, give it a nice thick stroke
													d3.select(this).style("stroke-width", '10px');
													hoverDiv.html("<b>" + timeOut + " survived:</b> " + d3.format(".1f")(x.invert(d3.mouse(this)[0])) + "<br/><b>probability of survival:</b> " + d3.format(".2%")(y.invert(d3.mouse(this)[1])));
													hoverDiv.transition()
															.duration(200)
															.style("opacity", 0.8)
															.style("left", (d3.mouse(this)[0]+90) + "px")
															.style("top", (d3.mouse(this)[1]+310) + "px");
												})
												.on("mousemove", function() {
													hoverDiv.html("<b>" + timeOut + " survived:</b> " + d3.format(".1f")(x.invert(d3.mouse(this)[0])) + "<br/><b>probability of survival:</b> " + d3.format(".2%")(y.invert(d3.mouse(this)[1])));
													hoverDiv.style("left", (d3.mouse(this)[0]+90) + "px")
															.style("top", (d3.mouse(this)[1]+310) + "px");
												})
												.on('mouseout', function () {
													d3.select(this).style("stroke-width", "3px");
													hoverDiv.transition().duration(500).style("opacity", 0);
												});
						lines.push(line);
					}
					
					updateLegend();
					
					updateRiskTable();
					
					updateSurvivalTables();
					
					defineSizes();
					
					sizeElements();
					
					calculatePositions();
					
					positionElements();
					
					makeLegendWhite();
					
					function updateLegend() {
						
						if(showLegend) {
							
							// add legend
							kaplanLegend = kaplanContainer.append("svg").attr("id", "kaplanLegend");
							
							// add labels and text to the legend
							kaplanLegend.selectAll('g')
										.data(survivalData)
										.enter()
										.append('g')
										.each(function(d, i) {
											var g = d3.select(this);
											g.append("rect")
												.attr("x", 0 + kaplanLegendPadding.left)
												.attr("y", i*kaplanLegendLineHeight + kaplanLegendPadding.top + kaplanLegendBoxMargin.top)
												.attr("class", "legend-box active")
												.attr("id", "legend-box-" + i)
												.attr("width", kaplanLegendBoxHeight)
												.attr("height", kaplanLegendBoxWidth)
												.style("fill", colors[i])
												.on("click", function() {
													// Determine if current line is visible
													var active = lines[i].active ? false : true;
													var newClass = active ? "inactive" : "active";
													// Hide or show the elements
													d3.select("#line-" + i).attr("class", newClass);
													d3.select("#legend-box-" + i).attr("class", "legend-box " + newClass);
													d3.select("#legend-text-" + i).attr("class", newClass);
													// Update whether or not the elements are active
													lines[i].active = active;
												})
												.on("mouseover", function() {
													d3.select("#line-" + i + ".active").style("stroke-width", '10px');
												})
												.on('mouseout', function () {
													d3.select("#line-" + i).style("stroke-width", '3px');
												});
											
											g.append("text")
												.attr("x", 0 + kaplanLegendPadding.left + kaplanLegendBoxWidth + kaplanLegendBoxMargin.right)
												.attr("y", i * kaplanLegendLineHeight + 10 + kaplanLegendPadding.top + kaplanLegendBoxMargin.top)
												.attr("class", "legend-text active")
												.attr("id", "legend-text-" + i)
												.attr("height", kaplanLegendLineHeight)
												.style("fill", colors[String(i)][1])
												.text(legendLabels[i])
											
										});
							
						}
						
					}
					
					function updateRiskTable() {
					
						if(showRiskTable) {
						
							var riskTableData = [];
							for(var i = 0; i < survivalData.length; i++) {
								riskTableData[i] = [];
								for(var j = 0; j < tickCountX; j++) {
									riskTableData[i][j] = 0;
								}
							}
							
							// lastIndex
							var lastIndex = [];
							for(var j = 0; j < survivalData.length; j++) {
								lastIndex[j] = 0;
							}
							
							for(var i = 0; i < tickCountX; i++) {
								
								var currentTimeTick = (maxXValue / tickCountX) * i;
								var nextTimeTick = (maxXValue / tickCountX) * (i + 1);
								if(nextTimeTick > maxXValue) {
									nextTimeTick = currentTimeTick;
								}
								
								for(var j = 0; j < survivalData.length; j++) {
									var survivalDataFound = false;
									// FIRST TICK
									if(i == 0) {
										riskTableData[j][i] = survivalData[j][0]["n"];
									} else {
										for(var k = lastIndex[j]; k < survivalData[j].length; k++) {
											if(k > 0) {
												lastIndex[j] = k-1;
											} else {
												k = 0;
											}
											if(survivalData[j][k]["t"] > currentTimeTick) {
												if(k > 0) {
													riskTableData[j][i] = survivalData[j][k-1]["n"] - survivalData[j][k-1]["d"];
												} else {
													riskTableData[j][i] = survivalData[j][k]["n"] - survivalData[j][k]["d"];
												}
												survivalDataFound = true;
												break;
											}
										}
										if(!survivalDataFound) {
											riskTableData[j][i] = survivalData[j][survivalData[j].length-1]["n"] - survivalData[j][survivalData[j].length-1]["d"];
										}
									}
								}
							}
							
							kaplanRiskTableContainer = kaplanContainer.append("foreignObject").attr("id", "kaplanRiskTableContainer");
							kapalnRiskTableContainerBody = kaplanRiskTableContainer.append("xhtml:body");
							kaplanRiskTable = kapalnRiskTableContainerBody.append("table").attr("id", "kaplanRiskTable").attr('style', 'border-collapse: collapse; table-layout: fixed; width: 100%;');
							for(var m = 0; m < riskTableData.length; m++) {
								var line = kaplanRiskTable.append('tr');
								if(m % 2 == 0) {
									line.attr('style', 'background-color: #D5D5D5;');
								} else {
									line.attr('style', 'background-color: #F2F2F2;');
								}
								for(var l = 0; l < tickCountX; l++) {
									line.append('td').html(riskTableData[m][l]).attr("style", "color: " + colors[m] + '; text-align: center; width: 4%;');
								}
							}
							
						}
						
					}
					
					function updateSurvivalTables() {
						
						// Update Table 1
						updateSurvivalTable1();
						
						// Update Table 2
						updateSurvivalTable2();
						
						function updateSurvivalTable1() {
							
							// Remove old table
							d3.select('#satistics1').remove();
							
							// Define data for new table
							var selected_subsets = scope.data.selected_subsets;
							var statistics1 = scope.data.statistics1;
							var header_statistics1 = [""];
							for(var i = 0; i < selected_subsets.length; i++) {
								header_statistics1.push(selected_subsets[i]);
							}
							
							if(mergeSubsets) {
								var merged = "";
								for(var j = 0; j < selectedSubsets.length; j++) {
									merged = merged + selectedSubsets[j] + "<br/>";
								}
								header_statistics1[1] = merged;
							}
							
							var heading = d3.select(root).append('h1').text("Cox Regression Result")
							var table = d3.select(root).append('table').attr('class', 'sr-survival-table').attr('id', 'statistics1');
							var thead = table.append('thead');
							var tbody = table.append('tbody');
							
							thead.append('tr')
								.selectAll('th')
								.data(header_statistics1)
								.enter()
								.append('th')
								.html(function(d) {
									return d;
								});
							
							for(var i = 0; i < statistics1.length; i++) {
								tbody.append('tr').attr('id', 'row-' + i);
								var row = tbody.select('#row-' + i);
								for(var property in statistics1[i]) {
									if(statistics1[i].hasOwnProperty(property)) {
										row.append('td').text(statistics1[i][property]);
									}
								}
							}
							
						}
						
						function updateSurvivalTable2() {
						
							// Remove old table
							d3.select('#satistics2').remove();
							
							// Define data for new table
							var header_statistics2 = ['Subset', 'Category', 'Cox Coefficient', 'Hazards Ratio', 'Standard Deviation', 'z-Score', 'Lower Range of Hazards Ratio <br> (95% Confidence Interval)', 'Upper Range of Hazards Ratio <br> (95% Confidence Interval)'];
							var statistics2 = scope.data.statistics2;
							
							var table = d3.select(root).append('table').attr('class', 'sr-survival-table').attr('id', 'statistics2');
							var thead = table.append('thead');
							var tbody = table.append('tbody');
							
							thead.append('tr')
								.selectAll('th')
								.data(header_statistics2)
								.enter()
								.append('th')
								.html(function(d) {
									return d;
								});
							
							for(var i = 0; i < statistics2.length; i++) {
								tbody.append('tr').attr('id', 'row-' + i);
								var row = tbody.select('#row-' + i);
								for(var property in statistics2[i]) {
									if(statistics2[i].hasOwnProperty(property)) {
										var value = statistics2[i][property];
										if(property == 0 && mergeSubsets) {
											var merged = "";
											for(var k = 0; k < selectedSubsets.length; k++) {
												merged = merged + selectedSubsets[k].replace("Subset ", "") + "<br/>";
											}
											value = merged;
										}
										if(property == 1) {
											if(value != "NO_CATEGORY_SELECTED") {
												if(value.includes("&")) {
													var splitted = value.trim().split(" & ");
													var merged = "";
													for(var j = 0; j < splitted.length; j++) {
														var cleaned_split = smartRUtils.shortenConcept(splitted[j]);
														if(merged == "") {
															merged = cleaned_split
														} else {
															merged = merged + " <br/> " + cleaned_split;
														}
													}
													value = merged;
												} else {
													value = smartRUtils.shortenConcept(new String(value));
												}
											}
										}
										row.append('td').html(value);
									}
								}
							}
							
						}
						
					}
					
					function defineSizes() {
						
						// Kaplan Plot
						kaplanPlotInnerHeight = plotHeight;
						kaplanPlotHeight = kaplanPlotInnerHeight + kaplanPlotPadding.top + kaplanPlotPadding.bottom + 2*kaplanPlotLabelHeight;
						kaplanPlotInnerWidth = plotWidth;
						kaplanPlotWidth = kaplanPlotInnerWidth + kaplanPlotPadding.left + kaplanPlotPadding.right;
						
						// Kaplan Legend
						if(showLegend) {
							kaplanLegendHeight = document.getElementById("kaplanLegend").getBBox().height + kaplanLegendPadding.top + kaplanLegendPadding.bottom + 5;
							kaplanLegendWidth = document.getElementById("kaplanLegend").getBBox().width + kaplanLegendPadding.left + kaplanLegendPadding.right;
						}
						
						// Kaplan Risk Table
						if(showRiskTable) {
							kaplanRiskTableWidth = (kaplanPlotInnerWidth / (tickCountX - 1)) * tickCountX;
							kaplanRiskTableHeight = kaplanLegendCount * 30;
							kaplanRiskTableContainerWidth = kaplanRiskTableWidth;
							kaplanRiskTableContainerHeight = kaplanRiskTableHeight + kaplanRiskTablePadding.top + kaplanRiskTablePadding.bottom;
						}
						
						// Kaplan Container
						// needs to be defined last because its size depends on the size of the other elements
						kaplanContainerWidth = kaplanPlotWidth + kaplanContainerPadding.left + kaplanContainerPadding.right;
						if(legendType == "outer" && showLegend) {
							kaplanContainerHeight = kaplanPlotHeight + kaplanLegendHeight + kaplanContainerPadding.top + kaplanContainerPadding.bottom;
						} else {
							kaplanContainerHeight = kaplanPlotHeight + kaplanContainerPadding.top + kaplanContainerPadding.bottom;
						}
						if(showRiskTable) {
							kaplanContainerHeight = kaplanContainerHeight + kaplanRiskTableContainerHeight;
						}
						
					}
					
					function sizeElements() {
						
						// Kaplan Container
						kaplanContainer.attr("width", kaplanContainerWidth).attr("height", kaplanContainerHeight);
						
						// Kaplan Plot
						kaplanPlot.attr("width", kaplanPlotWidth).attr("height", kaplanPlotHeight);
						
						// Kaplan Legend
						if(showLegend) {
							kaplanLegend.attr("width", kaplanLegendWidth).attr("height", kaplanLegendHeight);
						}
						
						// Kaplan Risk Table
						if(showRiskTable) {
							kaplanRiskTableContainer.attr("width", kaplanRiskTableContainerWidth).attr("height", kaplanRiskTableContainerHeight);
							kaplanRiskTable.attr("width", kaplanRiskTableWidth).attr("height", kaplanRiskTableHeight);
						}
						
					}
					
					function calculatePositions() {
						
						// Kaplan Container
						positions.kaplanContainerX = 0;
						positions.kaplanContainerY = 0;
						
						// Kaplan Plot
						positions.kaplanPlotX = kaplanContainerPadding.left;
						positions.kaplanPlotY = kaplanContainerPadding.top;
						if(legendType == "outer" && legendPosition == "top") {
							positions.kaplanPlotY = positions.kaplanPlotY + kaplanLegendHeight;
						}
						
						// Kaplan Plot X-Axis
						positions.kaplanPlotXAxisX = kaplanPlotPadding.left;
						positions.kaplanPlotXAxisY = kaplanPlotInnerHeight + kaplanPlotPadding.top;
						
						// Kaplan Plot X-Label
						positions.kaplanPlotXLabelX = kaplanPlotPadding.left + (kaplanPlotInnerWidth / 2);
						positions.kaplanPlotXLabelY = kaplanPlotHeight - kaplanPlotLabelHeight;
						
						// Kaplan Plot Y-Axis
						positions.kaplanPlotYAxisX = kaplanPlotPadding.left;
						positions.kaplanPlotYAxisY = kaplanPlotPadding.top;
						
						// Kaplan Plot Y-Label
						positions.kaplanPlotYLabelX = 0 - kaplanPlotPadding.top - (kaplanPlotInnerHeight / 2);
						positions.kaplanPlotYLabelY = kaplanPlotLabelHeight;
						
						// Kaplan Legend
						if(showLegend) {
							if(legendType == "inner") {
								if(legendPosition == "bottom") {
									positions.kaplanLegendX = kaplanContainerPadding.left + kaplanPlotPadding.left + kaplanLegendMargin.left;
									positions.kaplanLegendY = kaplanContainerPadding.top + kaplanPlotPadding.top + kaplanPlotInnerHeight - kaplanLegendHeight - kaplanLegendMargin.bottom;
								} else {
									positions.kaplanLegendX = kaplanContainerPadding.left + kaplanPlotWidth - kaplanLegendWidth - kaplanPlotPadding.right - kaplanLegendMargin.right;
									positions.kaplanLegendY = kaplanContainerPadding.top + kaplanPlotPadding.top + kaplanLegendMargin.top;
								}
							} else {
								positions.kaplanLegendX = kaplanContainerPadding.left + kaplanPlotPadding.left + ((kaplanPlotInnerWidth - kaplanLegendWidth)/2);
								if(legendPosition == "bottom") {
									positions.kaplanLegendY = kaplanPlotHeight + kaplanContainerPadding.bottom;
								} else {
									positions.kaplanLegendY = kaplanContainerPadding.top;
								}
							}
						}
						
						// Kaplan Plot Risk Table
						if(showRiskTable) {
							positions.kaplanRiskTableX = kaplanContainerPadding.left + kaplanPlotPadding.left - ((kaplanPlotInnerWidth/(tickCountX-1))/2);
							positions.kaplanRiskTableY = kaplanPlotHeight + kaplanContainerPadding.top + kaplanRiskTablePadding.top;
							if(showLegend && legendType == "outer") {
								positions.kaplanRiskTableY = positions.kaplanRiskTableY + kaplanLegendHeight;
							}
							positions.kaplanRiskTableContainerX = positions.kaplanRiskTableX;
							positions.kaplanRiskTableContainerY = positions.kaplanRiskTableY;
						}
						
					}
					
					function positionElements() {
						
						// Kaplan Container
						kaplanContainer.attr('x', positions.kaplanContainerX).attr('y', positions.kaplanContainerY);
						
						// Kaplan Plot (including axis & labels)
						kaplanPlot.attr('x', positions.kaplanPlotX).attr('y', positions.kaplanPlotY);
						theXAxis.attr('transform', 'translate(' + positions.kaplanPlotXAxisX + ',' + positions.kaplanPlotXAxisY + ')');
						theXLabel.attr('transform', 'translate(' + positions.kaplanPlotXLabelX + ',' + positions.kaplanPlotXLabelY + ')');
						theYAxis.attr('transform', 'translate(' + positions.kaplanPlotYAxisX + ',' + positions.kaplanPlotYAxisY + ')');
						theYLabel.attr('x', positions.kaplanPlotYLabelX).attr('y', positions.kaplanPlotYLabelY);
						
						// Kaplan Legend
						if(showLegend) {
							kaplanLegend.attr('x', positions.kaplanLegendX).attr('y', positions.kaplanLegendY);
						}
						
						// Kaplan Risk Table
						if(showRiskTable) {
							kaplanRiskTableContainer.attr('x', positions.kaplanRiskTableContainerX).attr('y', positions.kaplanRiskTableContainerY);
							kaplanRiskTable.attr('x', positions.kaplanRiskTableX).attr('y', positions.kaplanRiskTableY);
						}
						
					}
					
					// White Background for Legend
					function makeLegendWhite() {
						if(showLegend) {
							var strokeWidth = 2;
							kaplanLegend.insert("rect", ":first-child").attr("height", kaplanLegendHeight-2*strokeWidth).attr("width", kaplanLegendWidth-2*strokeWidth).attr("x", strokeWidth).attr("y", strokeWidth).attr("fill", "white").attr("stroke", "black").attr("stroke-width", strokeWidth);
						}
					}
					
					
				/* Drawing ends here */
				
			/* End d3.js */
			
		}
	
	}
	
]);
