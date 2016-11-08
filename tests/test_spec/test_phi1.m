function test_suite=test_phi1()
  initTestSuite;
end
function test_phi1_()
    S = jonswap; 
   S1=S; S1.S=S1.S.*phi1(S1.w,30);
end
