function test_suite=test_ftf()
  initTestSuite;
end
function test_ftf_()
    T = 1; 
   tp = dat2tp(load('sea.dat')); 
   RFC = tp2rfc(tp); 
   [t,F] = ftf(5.5e-10,cc2dam(RFC,5)/T,0.06,0.5);
end
