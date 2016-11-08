function test_suite=test_covplot()
  initTestSuite;
end
function test_covplot_()
    S1 = demospec;  
   S2 = demospec('dir'); 
   R1 = spec2cov(S1);  
   R2 = spec2cov(S2,0,30,2,100,[],2,[]); 
   subplot(211); L1 = 17; covplot(R1,L1,1,'.-'); 
   subplot(212); covplot(R2); 
 
   assert(R1.R(1:3)',... 
      [2.890364144652255, -0.547532218293498, -1.990124104043658], 1e-10); 
   close all;
end
