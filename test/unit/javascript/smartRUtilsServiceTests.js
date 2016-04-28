'use strict';

describe('smartRUtils', function() {

    var smartRUtils;

    beforeEach(module('smartRApp'));

    beforeEach(inject(function(_smartRUtils_) {
        smartRUtils = _smartRUtils_;
    }));

    it('must be defined', function() {
        expect(smartRUtils).toBeDefined();
    });

    it('has working conceptBoxMapToConceptKeys()', function() {
        var params = {a: {concepts: ['c1', 'c2'], valid: true}, 
            'foo bar __-!*()123 abc': {concepts: ['c3'], valid: true},
            b: {concepts: [], valid: true},
            c: {concepts: ['12--_- c31/??/*&^/foobar'], valid: true}
        };
        var expected = {a_n0: 'c1', a_n1: 'c2', 'foo bar __-!*()123 abc_n0': 'c3', c_n0: '12--_- c31/??/*&^/foobar'};
        var result = smartRUtils.conceptBoxMapToConceptKeys(params);
        expect(result).toEqual(expected);
    });

    it('has working makeSafeForCSS()', function() {
        var regex = /[_a-zA-Z]+[_a-zA-Z0-9-]*/;
        var testStr1 = '!@#$%^&*asdfghj _-;:,<.>/|\\"+=0987654321';
        var testStr2 = '__foo-bar-123';
        expect(testStr1.replace(regex, '')).not.toEqual('');
        expect(smartRUtils.makeSafeForCSS(testStr1).replace(regex, '')).toEqual('');
        expect(testStr2.replace(regex, '')).toEqual('');
        expect(smartRUtils.makeSafeForCSS(testStr2).replace(regex, '')).toEqual('');
    });

    it('has working shortenConcept()', function() {
        var testStr1 = '\\foo bar\\ab c d';
        var testStr2 = 'foo bar\\abcd';
        var testStr3 = '\\a\\b\\c\\d\\e\\f\\g\\h';
        var testStr4 = 'foo bar\\abcd\\';
        expect(smartRUtils.shortenConcept(testStr1)).toEqual('foo bar/ab c d');
        expect(smartRUtils.shortenConcept(testStr2)).toEqual('foo bar/abcd');
        expect(smartRUtils.shortenConcept(testStr3)).toEqual('g/h');
        expect(smartRUtils.shortenConcept(testStr4)).toEqual('foo bar/abcd');
    });
});
