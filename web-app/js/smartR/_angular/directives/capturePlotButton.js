//# sourceURL=capturePlotButton.js

'use strict';

window.smartRApp.directive('capturePlotButton', [function() {

    // aux for downloadSVG
    var copyWithCollapsedCSS = function(svgElement) {
        var relevantProperties = [
            'fill-opacity', 'fill', 'stroke', 'font-size', 'font-family',
            'shape-rendering', 'stroke-width'
        ];
        var clonedSvg = jQuery(svgElement).clone().attr('display', 'none');
        clonedSvg.insertAfter(svgElement);

        var cachedDefaults = {};
        var scratchSvg = jQuery(document.createElement('svg'))
            .attr('display', 'none')
            .appendTo(jQuery('body'));

        function getDefaultsForElement(jqElement) {
            var nodeName = jqElement.prop('nodeName');
            if (!cachedDefaults[nodeName]) {
                var newElement = jQuery(document.createElement(nodeName))
                    .appendTo(scratchSvg);

                cachedDefaults[nodeName] = window.getComputedStyle(newElement[0]);
            }
            return cachedDefaults[nodeName];
        }

        clonedSvg.find('*').each(function(idx, element) { // for each element in <svg>
            var computedStyle = window.getComputedStyle(element);

            var jqElem = jQuery(element);
            relevantProperties.forEach(function(property) { // for each property
                var effectiveStyle = computedStyle.getPropertyValue(property);
                var defaultStyle = getDefaultsForElement(jqElem).getPropertyValue(property);

                if (effectiveStyle !== defaultStyle) {
                    jqElem.attr(property, effectiveStyle);
                }
            });
        });

        scratchSvg.remove();

        return clonedSvg;
    };

    var downloadSVG = function(svgElement, fileName) {
        var serializer = new XMLSerializer();
        var clonedSvg = copyWithCollapsedCSS(svgElement);
        var xmlString = serializer.serializeToString(clonedSvg[0]);
        var blob = new Blob([xmlString], { type: 'image/svg+xml' });
        var svgBlobUrl = URL.createObjectURL(blob);
        var link = jQuery('<a/>')
            .attr('href', svgBlobUrl)
            .attr('download', fileName)
            .css('display', 'none');
        jQuery('body').append(link);
        link[0].click();
        link.remove();
        URL.revokeObjectURL(svgBlobUrl);
        clonedSvg.remove();
    };

    return {
        restrict: 'E',
        scope: {
            disabled: '=',
            filename: '@'
        },
        template:
            '<input type="button" value="Capture" class="heim-action-button" ng-click="capture()">',
        link: function(scope, elements) {

            var template_btn = elements.children()[0];
            template_btn.disabled = scope.disabled;

            scope.$watch('disabled', function (newValue) {
                template_btn.disabled = newValue;
            }, true);

            if (!scope.filename) {
                // default filename
                scope.filename = 'image.svg';
            }

            scope.capture = function() {
                var svgElement = jQuery('svg.visualization')[0];
                downloadSVG(svgElement, scope.filename);
            };

        }
    };
}]);
