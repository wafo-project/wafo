function test_suite=test_nt2lc()
  initTestSuite;
end
function test_nt2lc_()
    F0 = round(triu(rand(50),1)*10); 
   NT = cmat2nt(F0); 
   param = [-1 1 50]; 
   lc = nt2lc(param,NT); 
   plotlc(lc); 
 
   close all;
end
