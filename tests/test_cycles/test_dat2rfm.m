function test_suite=test_dat2rfm()
  initTestSuite;
end
function test_dat2rfm_()
    x = load('sea.dat'); 
   [RFM,u] = dat2rfm(x);    % Default parameters 
   subplot(1,2,1); cmatplot(u,u,RFM,3); 
   [RFM,u] = dat2rfm(x,0.5,[-2.5 2.5 50]); 
   subplot(1,2,2); cmatplot(u,u,RFM,3); 
 
   close all;
end
