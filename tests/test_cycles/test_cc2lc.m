function test_suite=test_cc2lc()
  initTestSuite;
end
function test_cc2lc_()
    tp = dat2tp(load('sea.dat')); 
   mM = tp2mm(tp); 
   lc = cc2lc(mM); 
 
   close all;
end
