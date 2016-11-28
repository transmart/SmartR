//# sourceURL=smartRUtils.js

'use strict';

window.smartRApp.factory('smartRUtils', ['$q', function($q) {

    var service = {};

    service.conceptBoxMapToConceptKeys = function smartRUtils_conceptBoxMapToConceptKeys(conceptBoxMap) {
        var allConcepts = {};
        Object.keys(conceptBoxMap).forEach(function(group) {
            var concepts = conceptBoxMap[group].concepts;
            concepts.forEach(function(concept, idx) {
                allConcepts[group + '_' + 'n' + idx] = concept;
            });
        });
        return allConcepts;
    };

    /**
     * Creates a CSS safe version of a given string
     * This should be used consistently across the whole of SmartR to avoid data induced CSS errors
     *
     * @param str
     * @returns {string}
     */
    service.makeSafeForCSS = function smartRUtils_makeSafeForCSS(str) {
        return String(str).replace(/[^a-z0-9]/g, function(s) {
            var c = s.charCodeAt(0);
            if (c === 32) {
                return '-';
            }
            if (c >= 65 && c <= 90) {
                return '_' + s.toLowerCase();
            }
            return '__' + ('000' + c.toString(16)).slice(-4);
        });
    };

    service.getElementWithoutEventListeners = function(cssSelector) {
        var element = document.getElementById(cssSelector);
        var copy = element.cloneNode(true);
        element.parentNode.replaceChild(copy, element);
        return copy;
    };

    service.shortenConcept = function smartRUtils_shortenConcept(concept) {
        var split = concept.split('\\');
        split = split.filter(function(str) { return str !== ''; });
        return split[split.length - 2] + '/' + split[split.length - 1];
    };

    $.fn.textWidth = function(text, css) {
        if (!$.fn.textWidth.fakeEl) {
            $.fn.textWidth.fakeEl = $('<span>').hide().appendTo(document.body);
        }
        var textEl = $.fn.textWidth.fakeEl.text(text || this.val() || this.text());
        if (css) {
            for (var key in css) {
                if (css.hasOwnProperty(key)) {
                    textEl.css(key, css[key]);
                }
            }
        } else {
            textEl.css('font', this.css('font'));
        }
        return $.fn.textWidth.fakeEl.width();
    };

    /**
     * This method creates a fake span to meassure the text with of a given string for given css attributes
     * @param text
     * @param font
     * @return text width in px
     */
    service.getTextWidth = function(text, font) {
        return $.fn.textWidth(text, font);
    };

    /**
     * This method returns the font-size for which the given text is still within boundaries of targetWidth
     * It also takes rotation into account
     *
     * @param text
     * @param css
     * @param startSize
     * @param targetWidth
     * @param rotation
     * @param shrinkStep
     * @return font size in px
     */
    service.scaleFont = function(text, css, startSize, targetWidth, rotation, shrinkStep)  {
        var _spaceGainedByRotation = function(r, a) {
           return r - r * Math.cos(a * Math.PI / 180);
        };

        var fontSize = startSize;
        while (fontSize > 0) { // should stop long before
            css['font-size'] = fontSize + 'px';
            var currentWidth = service.getTextWidth(text, css);
            if (currentWidth - _spaceGainedByRotation(currentWidth, rotation) < targetWidth) {
                break;
            }
            fontSize -= shrinkStep;
        }
        return fontSize;
    };

    /**
     * Executes callback with scroll position when SmartR mainframe is scrolled
     * @param function
     */
    service.callOnScroll = function(callback) {
        $('#sr-index').parent().scroll(function() {
            var scrollPos = $(this).scrollTop();
            callback(scrollPos);
        });
    };

    /**
     * Set minimum size for visualisation
     */
    service.prepareWindowSize = function(width, height) {
        $('#heim-tabs').css('min-width', parseInt(width) + 25);
        $('#heim-tabs').css('min-height', parseInt(height) + 25);
    };

    service.getScrollBarWidth = function() {
        var inner = document.createElement('p');
        inner.style.width = "100%";
        inner.style.height = "200px";

        var outer = document.createElement('div');
        outer.style.position = "absolute";
        outer.style.top = "0px";
        outer.style.left = "0px";
        outer.style.visibility = "hidden";
        outer.style.width = "200px";
        outer.style.height = "150px";
        outer.style.overflow = "hidden";
        outer.appendChild (inner);

        document.body.appendChild (outer);
        var w1 = inner.offsetWidth;
        outer.style.overflow = 'scroll';
        var w2 = inner.offsetWidth;
        if (w1 === w2) {
            w2 = outer.clientWidth;
        }

        document.body.removeChild (outer);

        return (w1 - w2);
    };

    /** 
     * removes all kind of elements that might live out of the viz directive (e.g. tooltips, contextmenu, ...)
     */
    service.cleanUp = function() {
        $('.d3-tip').remove();
    };

    service.countCohorts = function() {
        return !window.isSubsetEmpty(1) + !window.isSubsetEmpty(2);
    };

    service.getSubsetIds = function smartRUtil_getSubsetIds() {
        var defer = $q.defer();

        function resolveResult() {
            var res = window.GLOBAL.CurrentSubsetIDs.slice(1).map(function (v) {
                return v || null;
            });
            if (res.some(function (el) {
                return el !== null;
            })) {
                defer.resolve(res);
            } else {
                defer.reject();
            }
        }

        for (var i = 1; i <= window.GLOBAL.NumOfSubsets; i++) {
            if (!window.isSubsetEmpty(i) && !window.GLOBAL.CurrentSubsetIDs[i]) {
                window.runAllQueries(resolveResult, window.smartRPanel);
                return defer.promise;
            }
        }

        resolveResult();

        return defer.promise;
    };

    /**
     * Some cool Array functions
     */

    // fast unique()
    service.unique = function(arr, callback) {
        callback = typeof callback === 'undefined' ? function(d) { return d; } : callback;
        var a = [];
        var uniqArr = [];
        for (var i = 0, l = arr.length; i < l; i++) {
            var value = callback(arr[i]);
            if (a.indexOf(value) === -1) {
                a.push(value);
                uniqArr.push(arr[i]);
            }
        }
        return uniqArr;
    };

    service.getValuesForDimension = function(dimension, ascendingOrder) {
        var values = [];
        if (typeof ascendingOrder === 'undefined' || !ascendingOrder) {
            values = dimension.top(Infinity).map(function(record) { return dimension.accessor(record); });
        } else {
            values =dimension.bottom(Infinity).map(function(record) { return dimension.accessor(record); });
        }

        return values;
    };

    service.toggleLoadingScreen = function(visible) {
        document.getElementById('sr-loading').style.visibility = visible ? 'visible' : 'hidden';
    };

    service._intersectArrays = function(a, b) {
        var t;
        if (b.length > a.length) t = b, b = a, a = t; // indexOf to loop over shorter
        return a.filter(function (e) {
            if (b.indexOf(e) !== -1) return true;
        });
    };

    return service;
}]);
