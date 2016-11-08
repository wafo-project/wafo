function test_suite=test_cmatplot()
  initTestSuite;
end
function test_cmatplot_()
    param = [-1 1 64]; u=levels(param); 
   F = mktestmat(param,[-0.2 0.2],0.25,1/2); 
   cmatplot(F,1); 
   cmatplot(u,u,F,2); colorbar; 
   cmatplot(u,u,F,3); colorbar; 
   cmatplot(u,u,F,4); 
 
   close all;
end
