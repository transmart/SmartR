'use strict';

describe('fetchButton', function() {
    var $compile,
        $rootScope,
        $httpBackend,
        $q,
        element,
        smartRUtils,
        rServeService;

    beforeEach(module('smartRApp'));
    beforeEach(module('smartRTemplates'));

    beforeEach(inject(function(_$compile_, _$rootScope_, _$httpBackend_, _$q_, _smartRUtils_, _rServeService_) {
        $compile = _$compile_;
        $rootScope = _$rootScope_;
        $httpBackend = _$httpBackend_;
        $q = _$q_;
        smartRUtils = _smartRUtils_;
        rServeService = _rServeService_;

        $rootScope.conceptMap = {};

        element = $compile(angular.element('<fetch-button concept-map="conceptMap"></fetch-button>'))($rootScope);
        $rootScope.$digest();
    }));

    var _clickButton = function(loadDataIntoSessionReturn, executeSummaryStatsReturn) {
        var defer1 = $q.defer();
        var defer2 = $q.defer();
        spyOn(rServeService, 'loadDataIntoSession').and.returnValue(defer1.promise);
        spyOn(rServeService, 'executeSummaryStats').and.returnValue(defer2.promise);

        element.find('input').click();

        if (loadDataIntoSessionReturn) {
            eval('defer1.' + loadDataIntoSessionReturn);
            $rootScope.$digest();
        }


        if (executeSummaryStatsReturn) {
            eval('defer2.' + executeSummaryStatsReturn);
            $rootScope.$digest();
        }
    };

    var _prepareScope = function(cohorts, showSummaryStats, conceptMap) {
        spyOn(smartRUtils, 'countCohorts').and.returnValue(cohorts);
        element.isolateScope().showSummaryStats = showSummaryStats;
        element.isolateScope().conceptMap = conceptMap;
        $rootScope.$digest();
    };

    it('replaces element with content', function() {
        expect(element.find('input')).toBeDefined();
        expect(element.find('span')).toBeDefined();
    });

    it('should have the correct scope when clicked', function() {
        try { // we just want to see if the progress message is there
            element.find('input').click();
        } catch(e) {}
        expect(element.find('span').text()).toContain('Fetching data');
        expect(element.isolateScope().running).toBe(true);
        expect(element.isolateScope().loaded).toBe(false);
        expect(element.isolateScope().subsets).toEqual(0);
        expect(element.isolateScope().allSamples).toEqual(0);
    });

    it('should show another text after data is loaded if showSummaryStats is enabled', function() {
        _prepareScope(1, true, {foo: {concepts: ['concept'], valid: true}});
        try { // we just want to see if the progress message is there
            _clickButton('resolve()');
        } catch (e) {}
        expect(element.find('span').text()).toContain('Execute summary statistics');
        expect(element.isolateScope().running).toBe(true);
        expect(element.isolateScope().loaded).toBe(false);
        expect(element.isolateScope().subsets).toEqual(0);
        expect(element.isolateScope().allSamples).toEqual(0);
    });

    it('should have the correct scope when successful without showSummaryStats', function() {
        _prepareScope(2, false, {foo: {concepts: ['concept'], valid: true}});
        _clickButton('resolve()');
        expect(element.find('span').text()).toBe('Task complete! Go to the "Preprocess" or "Run Analysis" tab to continue.');
        expect(element.isolateScope().running).toBe(false);
        expect(element.isolateScope().loaded).toBe(true);
        expect(element.isolateScope().subsets).toEqual(2);
        expect(element.isolateScope().allSamples).toEqual(0);
    });

    it('should have the correct scope when successful with showSummaryStats', function() {
        _prepareScope(1, true, {foo: {concepts: ['concept'], valid: true}});
        _clickButton('resolve()', 'resolve({result: {allSamples: 1337, subsets: null}})');
        expect(element.find('span').text()).toBe('Task complete! Go to the "Preprocess" or "Run Analysis" tab to continue.');
        expect(element.isolateScope().running).toBe(false);
        expect(element.isolateScope().loaded).toBe(true);
        expect(element.isolateScope().subsets).toEqual(1);
        expect(element.isolateScope().allSamples).toEqual(1337);
    });

    it('should show error due to missing cohorts', function() {
        _prepareScope(0, false, {foo: {concepts: ['concept'], valid: true}});
        _clickButton();
        expect(element.find('span').text()).toBe('Error: No cohorts selected!');
    });

    it('should show error due to missing concepts', function() {
        _prepareScope(1, false, {foo: {concepts: [], valid: true}});
        _clickButton();
        expect(element.find('span').text()).toBe('Error: No concepts selected!');
    });
    
    it('should show error due to invalid concepts', function() {
        _prepareScope(1, false, {foo: {concepts: [], valid: false}});
        _clickButton();
        expect(element.find('span').text()).toBe('Error: Your data do not match the requirements! All fields must be green.');
    });
});

