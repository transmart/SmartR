'use strict';

// Avoid `console` errors in browsers that lack a console.
(function() {
    var method;
    var noop = function () {};
    var methods = [
        'assert', 'clear', 'count', 'debug', 'dir', 'dirxml', 'error',
        'exception', 'group', 'groupCollapsed', 'groupEnd', 'info', 'log',
        'markTimeline', 'profile', 'profileEnd', 'table', 'time', 'timeEnd',
        'timeline', 'timelineEnd', 'timeStamp', 'trace', 'warn'
    ];
    var length = methods.length;
    console = (window.console = window.console || {}); // jshint ignore:line

    while (length--) {
        method = methods[length];

        // Only stub undefined methods.
        if (!console[method]) {
            console[method] = noop;
        }
    }
}());

var smartRPanel = new Ext.Panel({
    id: 'smartRPanel',
    title: 'SmartR',
    region: 'center',
    split: true,
    height: 90,
    layout: 'fit',
    collapsible: true,
    autoScroll: true,
    tbar: new Ext.Toolbar({
        id: 'smartRToolbar',
        title: 'R Scripts',
        items: []
    }),
    listeners: {
        render: function (panel) {
            /**
        * WORKAROUND : code below is needed to reorder the javascript script load that're broken due to
        * ExtJS panel
        */
            // start workaround
            var updater = panel.getUpdater();
            updater.on('update', function() {
                var panelBody = jQuery(arguments[0].dom);
                var scripts = panelBody.children('script');
                scripts.remove(); // remove scripts from panel body
                panelBody.append(scripts); // re-append again
            });
            updater.update({
                url: pageInfo.basePath + '/smartR/index',
                method: 'POST',
                scripts: false });
            // end workaround
        }
    }
});

window.addSmartRPanel = function addSmartRPanel(parentPanel) {
    parentPanel.insert(4, smartRPanel);
};

function cleanUpSmartR() {
    var d3tips = document.getElementsByClassName('d3-tip');
    Array.prototype.forEach.call(d3tips, function(el) { // d3tips is array-like object
        el.parentNode.removeChild(el);
    });
}
cleanUpSmartR();

