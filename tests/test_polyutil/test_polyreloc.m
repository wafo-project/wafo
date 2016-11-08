function test_suite=test_polyreloc()
  initTestSuite;
end
function test_polyreloc_()
    assert(polyreloc([1,1], 2), [1, -1], eps); 
   assert(polyreloc([1,1], 4), [1, -3], eps); 
   assert(polyreloc([1,1], 2, 2), [1,1], eps);
end
