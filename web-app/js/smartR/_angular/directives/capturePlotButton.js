//# sourceURL=capturePlotButton.js

'use strict';

window.smartRApp.directive('capturePlotButton', [function() {

    return {
        restrict: 'E',
        scope: {
            disabled: '=?',
            filename: '@',
            target: '@'
        },
        template: '<button id="sr-capture-button" ng-disabled="disabled">Capture SVG</button>',
        link: function(scope, elements) {
            if (!scope.filename) {
                // default filename
                scope.filename = 'image.svg';
            }


            var template_btn = elements.children()[0];
            template_btn.addEventListener('click', function() {
                var svg = $(scope.target + ' svg');
                if (!svg.length) {
                    return;
                }

                var smartRSheets = [];
                for (var i = 0; i < document.styleSheets.length; i++) {
                    var sheet = document.styleSheets[i];
                    if (sheet.href && sheet.href.indexOf('smart-r') !== -1) {
                        smartRSheets.push(sheet);
                    }
                }

                var rules = [];
                smartRSheets.forEach(function(d) {
                    for (var key in d.cssRules) {
                        if (d.cssRules.hasOwnProperty(key) && d.cssRules[key] instanceof CSSStyleRule) {
                            rules.push(d.cssRules[key].cssText);
                        }
                    }
                });

                function b64toBlob(b64Data, contentType, sliceSize) {
                    contentType = contentType || '';
                    sliceSize = sliceSize || 512;

                    var byteCharacters = atob(b64Data);
                    var byteArrays = [];

                    for (var offset = 0; offset < byteCharacters.length; offset += sliceSize) {
                        var slice = byteCharacters.slice(offset, offset + sliceSize);

                        var byteNumbers = new Array(slice.length);
                        for (var i = 0; i < slice.length; i++) {
                            byteNumbers[i] = slice.charCodeAt(i);
                        }

                        var byteArray = new Uint8Array(byteNumbers);

                        byteArrays.push(byteArray);
                    }

                    var blob = new Blob(byteArrays, {type: contentType});
                    return blob;
                }

                var defs = '<defs><style type="text/css"><![CDATA[' + rules.join('') + ']]></style></defs>';
                svg.attr({version: '1.1' , xmlns:"http://www.w3.org/2000/svg"});
                svg.append(defs);
                $(scope.target + ' svg').wrap('<div id="sr-capture-container"></div>');
                var b64 = btoa(unescape(encodeURIComponent($('#sr-capture-container').html())));
                var blob = b64toBlob(b64, 'image/svg+xml');
                var blobUrl = URL.createObjectURL(blob);
                var a = document.createElement('a');
                a.style = 'display: none';
                a.href = blobUrl;
                a.download = scope.filename;
                a.click();
                URL.revokeObjectURL(blobUrl);
                $(scope.target + ' svg').unwrap();
                $(scope.target + ' svg defs').remove();
            });
        }
    };
}]);
