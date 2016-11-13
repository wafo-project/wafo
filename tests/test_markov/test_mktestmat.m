function test_suite=test_mktestmat()
  initTestSuite;
end
function test_mktestmat_()
    [F,Fh] = mktestmat([-1 1 32],[-0.2 0.2], 0.25,1/2); 
   u = levels([-1 1 32]);  
   cmatplot(u,u,F,3); axis('square'); 
   [F,Fh] = mktestmat([-1 1 32],[-0.2 0.2], 0.25,1/2,-Inf); 
   cmatplot(u,u,F,3); axis('square'); 
 
   close all;
end
