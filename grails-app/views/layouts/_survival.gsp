
<script type="text/ng-template" id="survival">
	
	<div ng-controller="SurvivalController">
		
		<tab-container>
			
			<workflow-tab tab-name="Fetch Data" disabled="fetch.disabled">
				<concept-box style="display: inline-block;"
							 concept-group="fetch.conceptBoxes.time"
							 type="LD-numerical"
							 min="1"
							 max="1"
							 label="Time"
							 tooltip="Select time variable from the Data Set Explorer Tree and drag it into the box. For example, 'Survival Time'. This variable is required.">
				</concept-box>
				<concept-box style="display: inline-block;"
							 concept-group="fetch.conceptBoxes.category"
							 type="LD-categorical"
							 id="categories"
							 min="0"
							 max="4"
							 label="Category (optional)"
							 tooltip="Select a variable on which you would like to group the cohort and drag it into the box. For example, 'Cancer Stage'.">
				</concept-box>
				<concept-box style="display: inline-block;"
							 concept-group="fetch.conceptBoxes.censoring"
							 type="LD-categorical"
							 id="censoring"
							 min="0"
							 max="1"
							 label="Censoring Variable (optional)"
							 tooltip="Drag the item for which to perform censoring in the analysis into this box. For example, when performing Overall survival analysis, drag 'Survival status = alive' into this box. This variable is not obligatory. ">
				</concept-box>
				<br/>
				<fetch-button concept-map="fetch.conceptBoxes"
							  loaded="fetch.loaded"
							  running="fetch.running"
							  allowed-cohorts="[1, 2]">
				</fetch-button>
			</workflow-tab>
			
			<workflow-tab tab-name="Run Analysis" disabled="runAnalysis.disabled">
				<div class="display-aside">
					<div class="heim-input-field sr-input-area">
						<h2>Plot Width <i class="ui-icon ui-icon-info sr-tooltip-dialog" title="Defines the Width of the survival plot."></i></h2>
						<fieldset class="heim-radiogroup">
							<label>
								<input type="number" ng-model="runAnalysis.params.plotWidth" value="800" min="350" max="900">
							</label>
						</fieldset>
					</div>
					<div class="heim-input-field sr-input-area">
						<h2>Plot Height <i class="ui-icon ui-icon-info sr-tooltip-dialog" title="Defines the Height of the survival plot."></i></h2>
						<fieldset class="heim-radiogroup">
							<label>
								<input type="number" ng-model="runAnalysis.params.plotHeight" value="500" min="350" max="600">
							</label>
						</fieldset>
					</div>
					<div class="heim-input-field sr-input-area">
						<h2>Time in <i class="ui-icon ui-icon-info sr-tooltip-dialog" title="Defines the domain the data is given in."></i></h2>
						<fieldset class="heim-radiogroup">
							<label>
								<input type="radio" ng-model="runAnalysis.params.timeIn" value="days" checked> Days
							</label>
							<label>
								<input type="radio" ng-model="runAnalysis.params.timeIn" value="months"> Months
							</label>
							<label>
								<input type="radio" ng-model="runAnalysis.params.timeIn" value="years"> Years
							</label>
						</fieldset>
					</div>
					<div class="heim-input-field sr-input-area">
						<h2>Time out <i class="ui-icon ui-icon-info sr-tooltip-dialog" title="Defines the domain the data should be displayed in the graph."></i></h2>
						<fieldset class="heim-radiogroup">
							<label>
								<input type="radio" ng-model="runAnalysis.params.timeOut" value="days" checked> Days
							</label>
							<label>
								<input type="radio" ng-model="runAnalysis.params.timeOut" value="months"> Months
							</label>
							<label>
								<input type="radio" ng-model="runAnalysis.params.timeOut" value="years"> Years
							</label>
						</fieldset>
					</div>
				</div>
				<div class="display-aside">
					<div class="heim-input-field sr-input-area">
						<h2>Legend Type <i class="ui-icon ui-icon-info sr-tooltip-dialog" title="Defines if the legend will be displayed inside or outside of the plot."></i></h2>
						<fieldset class="heim-radiogroup">
							<label>
								<input type="radio" ng-model="runAnalysis.params.legendType" value="inner" checked> Inner
							</label>
							<label>
								<input type="radio" ng-model="runAnalysis.params.legendType" value="outer"> Outer
							</label>
						</fieldset>
					</div>
					<div class="heim-input-field sr-input-area">
						<h2>Legend Position <i class="ui-icon ui-icon-info sr-tooltip-dialog" title="Defines the location the legend will be displayed."></i></h2>
						<fieldset class="heim-radiogroup">
							<label>
								<input type="radio" ng-model="runAnalysis.params.legendPosition" value="bottom" checked> Bottom
							</label>
							<label>
								<input type="radio" ng-model="runAnalysis.params.legendPosition" value="top"> Top
							</label>
						</fieldset>
					</div>
					<div class="heim-input-field sr-input-area">
						<h2>Show Risk Table <i class="ui-icon ui-icon-info sr-tooltip-dialog" title="Defines if the risk table will be displayed."></i></h2>
						<fieldset class="heim-radiogroup">
							<label>
								<input type="radio" ng-model="runAnalysis.params.showRiskTable" value="TRUE" checked> True
							</label>
							<label>
								<input type="radio" ng-model="runAnalysis.params.showRiskTable" value="FALSE"> False
							</label>
						</fieldset>
					</div>
				</div>
				<div class="display-aside">
					<div class="heim-input-field sr-input-area">
						<h2>Merge Subsets <i class="ui-icon ui-icon-info sr-tooltip-dialog" title="Merges the selected subsets by the selected categories. Results in one graph per category."></i></h2>
						<fieldset class="heim-radiogroup" ng-disabled="common.subsets < 2">
							<label>
								<input type="radio" ng-model="runAnalysis.params.mergeSubsets" value="FALSE" checked> False
							</label>
							<label>
								<input type="radio" ng-model="runAnalysis.params.mergeSubsets" value="TRUE"> True
							</label>
						</fieldset>
					</div>
					<div class="heim-input-field sr-input-area">
						<h2>Merge Categories  <i class="ui-icon ui-icon-info sr-tooltip-dialog" title="Merges the selected categories. Results in one graph per subset."></i></h2>
						<fieldset class="heim-radiogroup" ng-disabled="common.categories < 2">
							<label>
								<input type="radio" ng-model="runAnalysis.params.mergeCategories" value="FALSE" checked> False
							</label>
							<label>
								<input type="radio" ng-model="runAnalysis.params.mergeCategories" value="TRUE"> True
							</label>
						</fieldset>
					</div>
					<div class="heim-input-field sr-input-area">
						<h2>Censor Interpretation  <i class="ui-icon ui-icon-info sr-tooltip-dialog" title="Defines the interpretation of the specified censoring event."></i></h2>
						<fieldset class="heim-radiogroup" ng-disabled="common.censoring < 1">
							<label>
								<input type="radio" ng-model="runAnalysis.params.censorInterpretation" value="negative" checked> Event is 1
							</label>
							<label>
								<input type="radio" ng-model="runAnalysis.params.censorInterpretation" value="positive"> Event is 0
							</label>
						</fieldset>
					</div>
				</div>
				<hr class="sr-divider">
				<run-button button-name="Create Plot"
							store-results-in="runAnalysis.scriptResults"
							script-to-run="run"
							arguments-to-use="runAnalysis.params"
							running="runAnalysis.running">
				</run-button>
				<capture-plot-button filename="survivalplot.svg" disabled="runAnalysis.download.disabled" target="survival-plot"></capture-plot-button>
				<br/>
				<br/>
				<survival-plot data="runAnalysis.scriptResults" width="1000" height="700"></survival-plot>
			</workflow-tab>
			
		</tab-container>
		
	</div>
	
</script>
