function test_suite=test_tmaspec()
  initTestSuite;
end
function test_tmaspec_()
   S = tmaspec(3,[0 0 1], inf);   % Bretschneider spectrum Hm0=7, Tp=11 
  S2 = bretschneider(3); 
  assert(S.S, S2.S, 1e-10)  
  S = tmaspec(1.5,[0 0 1], inf);  % The same, cut at wc = 1.5 
  S2 = bretschneider(1.5); 
  assert(S.S, S2.S, 1e-10)
end
