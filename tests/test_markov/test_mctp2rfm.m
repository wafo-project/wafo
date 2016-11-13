function test_suite=test_mctp2rfm()
  initTestSuite;
end
function test_mctp2rfm_()
    param = [-1 1 32]; u = levels(param); 
   F = mktestmat(param,[-0.2 0.2],0.15,2); 
   Frfc = mctp2rfm({F []}); 
   cmatplot(u,u,Frfc); 
 
   close all;
end
