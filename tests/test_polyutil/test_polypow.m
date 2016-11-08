function test_suite=test_polypow()
  initTestSuite;
end
function test_polypow_()
    assert(polypow([1,1], 2), [1,2,1], 1e-12); 
   assert(polypow([1,1], 3), [1,3,3,1], 1e-12);
end
