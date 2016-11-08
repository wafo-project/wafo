function test_suite=test_tp2arfc()
  initTestSuite;
end
function test_tp2arfc_()
    x = load('sea.dat'); tp=dat2tp(x); 
   ARFC=tp2arfc(tp);      % Default (min-to-Max cycles in residual) 
   ccplot(ARFC); 
 
   close all;
end
