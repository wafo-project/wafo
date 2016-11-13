function test_suite=test_smc2rfm()
  initTestSuite;
end
function test_smc2rfm_()
    P = [0.9, 0.1;0.05 0.95]; 
   param = [-1 1 32]; u = levels(param); 
   [F1,F2] = mktestmat(param,[-0.3 0.3],0.15,1,-Inf); 
   Frfc = smc2rfm(P,{F1;F2}); 
   cmatplot(u,u,Frfc) 
 
   close all;
end
