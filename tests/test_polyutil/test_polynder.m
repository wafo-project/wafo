function test_suite=test_polynder()
  initTestSuite;
end
function test_polynder_()
    assert(polynder([3, 2, 1],-1), [1,1,1,0], 1e-10); 
   assert(polynder([1,1,1],0), [1,1,1], 1e-10); 
   assert(polynder([3,2,1],1), [6, 2], 1e-10); 
   assert(polynder([12,6,2],2), 24, 1e-10);
end
