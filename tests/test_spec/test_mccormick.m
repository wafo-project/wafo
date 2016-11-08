function test_suite=test_mccormick()
  initTestSuite;
end
function test_mccormick_()
    S = mccormick(3,[6.5 10]); 
   plotspec(S) 
   assert( S.S(60:63), [7.30656984162512,6.96602463131717,... 
                        6.61130359201946, 6.25051761036712]', 1e-10)
end
