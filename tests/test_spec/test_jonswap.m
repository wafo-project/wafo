function test_suite=test_jonswap()
  initTestSuite;
end
function test_jonswap_()
    S = jonswap(3,[0 0 1]); 
   S2 = bretschneider(3); 
   assert(S.S, S2.S, 1e-10)
end
