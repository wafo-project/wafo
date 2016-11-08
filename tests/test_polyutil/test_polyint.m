function test_suite=test_polyint()
  initTestSuite;
end
function test_polyint_()
   assert(polyint([3,2,1]), [1,1,1,0], 1e-10); 
  assert(polyint([3,2,1], 1), 3, 1e-10); 
  assert(polyint([3,2,1], -1, 1), 4, 1e-10);
end
