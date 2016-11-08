function test_suite=test_specnorm()
  initTestSuite;
end
function test_specnorm_()
    S = jonswap; 
   [Sn,mn4] = specnorm(S); 
   assert(spec2mom(Sn,2), [1,1], 1e-4)     % Should be equal to one!
end
