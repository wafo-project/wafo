function test_suite=test_dtp2rfm()
  initTestSuite;
end
function test_dtp2rfm_()
    x = load('sea.dat');                   % Load data 
   [dtp,u,tp] = dat2dtp([-2 2 32],x,0.2); % Discrete TP & rainflow filter 0.2 
   RFM = dtp2rfm(dtp,32);                 % Calculate rainflow matrix 
   cmatplot(u,u,RFM,3); colorbar;         % Plot rainflow matrix 
 
   close all;
end
