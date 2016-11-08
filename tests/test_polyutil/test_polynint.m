function test_suite=test_polynint()
  initTestSuite;
end
function test_polynint_()
    assert(polynint([1/2, 1, 1],-1), [1,1], 1e-10); 
   assert(polynint([1,1,1],0), [1,1,1], 1e-10); 
   assert(polynint([3,2,1],1), [1,1,1,0], 1e-10); 
   assert(polynint([12,6,2],2), [1,1,1,0,0], 1e-10);
end
