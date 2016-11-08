function test_suite=test_dtp2arfm()
  initTestSuite;
end
function test_dtp2arfm_()
    x = load('sea.dat');                   % Load data 
   [dtp,u,tp] = dat2dtp([-2 2 32],x,0.2); % Discrete TP & rainflow filter 0.2 
   ARFM = dtp2arfm(dtp,32);            % Calculate asymmetric rainflow matrix 
   cmatplot(u,u,ARFM,3); colorbar;      % Plot rainflow matrix 
 
   close all;
end
