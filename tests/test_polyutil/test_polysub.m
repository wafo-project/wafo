function test_suite=test_polysub()
  initTestSuite;
end
function test_polysub_()
    assert(polysub([1,2,3],[1,2]), [1,1,1], 1e-10);
end
