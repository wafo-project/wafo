function test_suite=test_plotsn()
  initTestSuite;
end
function test_plotsn_()
    sn = load('sn.dat'); s = sn(:,1); N = sn(:,2); 
   [e,beta,s2] = plotsn(s,N,2);   % S-N, x-axis = S, y-axis = N 
   [e,beta,s2] = plotsn(s,N,12);  % N-S, x-axis = N, y-axis = S 
 
   close all;
end
