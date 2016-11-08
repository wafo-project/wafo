function test_suite=test_res2arfc()
  initTestSuite;
end
function test_res2arfc_()
    x = load('sea.dat'); tp=dat2tp(x); 
   [ARFC,res]=tp2arfc4p(tp);      % Default (min-to-Max cycles in residual) 
   ARFC_res = res2arfc(res);      % Cycles in residual 
   plotcc(ARFC); hold on; plot(ARFC_res(:,1),ARFC_res(:,2),'r.'); hold off; 
 
   close all;
end
