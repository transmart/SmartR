//# sourceURL=conceptBox.js

'use strict';

window.smartRApp.directive('conceptBox', [
    '$rootScope',
    '$http',
    function($rootScope, $http) {
        return {
            restrict: 'E',
            scope: {
                conceptGroup: '=',
                label: '@',
                tooltip: '@',
                min: '@',
                max: '@',
                type: '@'
            },
            templateUrl: $rootScope.smartRPath +  '/js/smartR/_angular/templates/conceptBox.html',
            link: function(scope, element) {
                var max = parseInt(scope.max);
                var min = parseInt(scope.min);

                var template_box = element[0].querySelector('.sr-drop-input'),
                    template_btn = element[0].querySelector('.sr-drop-btn'),
                    template_tooltip = element[0].querySelector('.sr-tooltip-dialog');

                // instantiate tooltips
                $(template_tooltip).tooltip({track: true, tooltipClass:"sr-ui-tooltip"});

                var _clearWindow = function() {
                    $(template_box).children().remove();
                };

                var _getConcepts = function() {
                    return $(template_box).children().toArray().map(function(childNode) {
                        return childNode.getAttribute('conceptid');
                    });
                };

                var _activateDragAndDrop = function() {
                    var extObj = Ext.get(template_box);
                    var dtgI = new Ext.dd.DropTarget(extObj, {ddGroup: 'makeQuery'});
                    dtgI.notifyDrop = dropOntoCategorySelection; // jshint ignore:line
                };

                var typeMap = {
                    hleaficon: 'HD',
                    null: 'LD-categorical', // FIXME: alphaicon does not exist yet in transmartApp master branch
                    valueicon: 'LD-numerical'
                };
                var _containsOnlyCorrectType = function() {
                    return $(template_box).children().toArray().every(function(childNode) {
                        return typeMap[childNode.getAttribute('setnodetype')] === scope.type;
                    });
                };

                var _getNodeDetails = function(conceptKeys, callback) {
                    var request = $http({
                        url: pageInfo.basePath + '/SmartR/nodeDetails',
                        method: 'POST',
                        config: {
                            timeout: 10000
                        },
                        data: {
                            conceptKeys: conceptKeys
                        }
                    });

                    request.then(
                        callback,
                        function() {
                            alert('Could not fetch node details. Network connection lost?');
                        });
                };

                // activate drag & drop for our conceptBox and color it once it is rendered
                scope.$evalAsync(function() {
                    _activateDragAndDrop();
                });

                // bind the button to its clearing functionality
                template_btn.addEventListener('click', function() {
                    _clearWindow();
                });

                // this watches the childNodes of the conceptBox and updates the model on change
                new MutationObserver(function() {
                    scope.conceptGroup.concepts = _getConcepts(); // update the model
                    scope.validate();
                    scope.$apply();
                }).observe(template_box, { childList: true });

                scope.validate = function() {
                    scope.instructionMinNodes = scope.conceptGroup.concepts.length < min;
                    scope.instructionMaxNodes = max !== -1 && scope.conceptGroup.concepts.length > max;
                    scope.instructionNodeType = !_containsOnlyCorrectType();
                    if (scope.type === 'HD' && scope.conceptGroup.concepts.length > 1) {
                        _getNodeDetails(scope.conceptGroup.concepts, function(response) {
                            if (Object.keys(response.data).length < 2) {
                                var platforms = response.data[Object.keys(response.data)[0]].platforms;
                                scope.instructionNodePlatform = !platforms.every(function(el) { 
                                    return el.title === platforms[0].title;
                                });
                            } else {
                                scope.instructionNodePlatform = true;
                            }
                        });
                    } else {
                        scope.instructionNodePlatform = false;
                    }
                };

                scope.$watchGroup([
                    'instructionNodeType',
                    'instructionNodePlatform',
                    'instructionMaxNodes',
                    'instructionMinNodes'],
                    function(newValue) {
                        var instructionNodeType = newValue[0],
                            instructionNodePlatform = newValue[1],
                            instructionMaxNodes = newValue[2],
                            instructionMinNodes = newValue[3];
                        scope.conceptGroup.valid = !(instructionNodeType ||
                                                     instructionNodePlatform ||
                                                     instructionMaxNodes ||
                                                     instructionMinNodes);
                    });

                scope.validate();
            }
        };
    }]);
