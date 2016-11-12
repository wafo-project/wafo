function test_suite=test_cmat2dmat()
  initTestSuite;
end
function test_cmat2dmat_()
    param = [-1 1 32];  
   F = mktestmat(param); 
   Dmat = cmat2dmat(param,F,6); 
   cmatplot(Dmat); 
 
   close all;
end
