function test_suite=test_genchol()
  initTestSuite;
end
function test_genchol_()
   H = hilb(10); 
  tol   = 1e-6; 
  [L,P] = genchol(H,tol); 
  spy(L*L.'-H(P,P)); 
 
 close all;
end
