function test_suite=test_ratio()
  initTestSuite;
end
function test_ratio_()
   assert(ratio(2,1,1,1), cosh(2)/cosh(1), 1e-10) 
  assert(ratio(2,1,1,-1), cosh(2)/sinh(1), 1e-10) 
  assert(ratio(2,1,-1,1), sinh(2)/cosh(1), 1e-10) 
  assert(ratio(2,1,-1,-1), sinh(2)/sinh(1), 1e-10)
end
