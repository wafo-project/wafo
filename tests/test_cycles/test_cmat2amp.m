function test_suite=test_cmat2amp()
  initTestSuite;
end
function test_cmat2amp_()
    x = load('sea.dat');                   % Load data 
   [dtp,u,tp] = dat2dtp([-2 2 32],x,0.2); % Discrete TP & rainflow filter 0.2 
   RFM = dtp2rfm(dtp,32);                 % Calculate rainflow matrix 
   amp_hist = cmat2amp([-2 2 32],RFM);    % Get amplitude histigram 
   bar(amp_hist(:,1),amp_hist(:,2));      % Plot histogram 
 
   close all;
end
