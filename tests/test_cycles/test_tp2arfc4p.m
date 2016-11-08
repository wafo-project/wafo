function test_suite=test_tp2arfc4p()
  initTestSuite;
end
function test_tp2arfc4p_()
    x = load('sea.dat');  
   tp = dat2tp(x); 
   [ARFC, res]=tp2arfc4p(tp);  % Default (min-to-Max cycles in residual) 
   ccplot(ARFC); 
   % res 
 
   close all;
end
