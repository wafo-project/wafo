function test_suite=test_polyrescl()
  initTestSuite;
end
function test_polyrescl_()
    assert(polyrescl([1,1], 2), [0.5, 1], eps); 
   assert(polyrescl([1,1], 4), [0.25, 1], eps); 
   assert(polyrescl([1,1], 2, 2), [1, 2], eps);
end
