function test_suite=test_polyadd()
  initTestSuite;
end
function test_polyadd_()
   assert(polyadd([1,2,3],[1,2]), [1,3,5], 1e-10);
end
