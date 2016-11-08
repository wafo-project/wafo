function test_suite=test_normspec()
  initTestSuite;
end
function test_normspec_()
    S = jonswap; 
   [Sn,mn4] = normspec(S); 
   assert(spec2mom(Sn,2), [1,1], 1e-4)     % Should be equal to one!
end
