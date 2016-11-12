function test_suite=test_cmat2dam()
  initTestSuite;
end
function test_cmat2dam_()
    param = [-1 1 32];  
   F = mktestmat(param); 
   bv = 3:8;  
   D = cmat2dam(param,F,bv);  
   plot(bv,D,'x-'); 
 
   close all;
end
