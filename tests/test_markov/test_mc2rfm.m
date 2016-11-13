function test_suite=test_mc2rfm()
  initTestSuite;
end
function test_mc2rfm_()
    F = magic(3); 
   Q = mat2tmat(F); 
   Frfc = mc2rfm(Q); 
   assert(Frfc, [0.0,   0.006666666666667, 0.148888888888889;... 
                 0.0,   0.0,   0.140;... 
                 0.0,   0.0,   0.0], 1e-10)
end
