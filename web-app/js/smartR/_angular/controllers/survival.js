//# sourceURL=survival.js

'use strict';

window.smartRApp.controller('SurvivalController', [
	'$scope',
	'smartRUtils',
	'commonWorkflowService',
	function($scope, smartRUtils, commonWorkflowService) {
		
		commonWorkflowService.initializeWorkflow('survival', $scope);
		
		$scope.fetch = {
			disabled: false,
			running: false,
			loaded: false,
			conceptBoxes: {
				time: {concepts: [], valid: false},
				category: {concepts: [], valid: true},
				censoring: {concepts: [], valid: true}
			}
		};
		
		$scope.common = {
			subsets: 0,
			categories: 0
		};
		
		$scope.runAnalysis = {
			disabled: true,
			running: false,
			params: {
				plotWidth: '800',
				plotHeight: '500',
				timeIn: 'days',
				timeOut: 'days',
				legendType: 'inner',
				legendPosition: 'bottom',
				showRiskTable: 'TRUE',
				mergeSubsets: 'FALSE',
				mergeCategories: 'FALSE',
				censorInterpretation: 'negative'
			},
			download: {
				disabled: true,
			},
			scriptResults: {}
		};
		
		$scope.$watchGroup(['fetch.running', 'runAnalysis.running'],
			function(newValues) {
				var fetchRunning = newValues[0],
					runAnalysisRunning = newValues[1];
				
				// clear old results
				if (fetchRunning) {
					$scope.runAnalysis.scriptResults = {};
					$scope.common.subsets = smartRUtils.countCohorts();
					$scope.common.categories = countCategories();
					$scope.common.censoring = countCensoring();
				}
				
				// disable tabs when certain criteria are not met
				$scope.fetch.disabled = runAnalysisRunning;
				$scope.runAnalysis.disabled = fetchRunning || !$scope.fetch.loaded;
				
				// disable buttons when certain criteria are not met
				$scope.runAnalysis.download.disabled = runAnalysisRunning || $.isEmptyObject($scope.runAnalysis.scriptResults);
			}
		);
		
	}]);
	
	
function countCategories() {
	return document.querySelectorAll('#categories .conceptUnselected').length;
}

function countCensoring() {
	return document.querySelectorAll('#censoring .conceptUnselected').length;
}