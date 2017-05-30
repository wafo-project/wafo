function test_suite=test_finddep()
  initTestSuite;
end
function test_finddep_()
     [C,W,M] = finddep(fullfile(waforoot, 'fileutil')); 
    assert(W.m(1:4), {'Contents.m','bindiff.m','cdtomfile.m','diary2m.m'}); 
    % assert(size(C{1,1}), [52,52]); 
    assert(C{1,1}(1:4,1:4), sparse([2,3,4],[2,3,4], ones(1, 3)));
end
