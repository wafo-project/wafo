function test_suite=test_spec2cov2()
  initTestSuite;
end
function test_spec2cov2_()
       S = jonswap; 
      dt = 0.1;   
      R = spec2cov2(S,3,256,dt);
end
