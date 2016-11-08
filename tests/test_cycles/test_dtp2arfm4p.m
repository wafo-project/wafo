function test_suite=test_dtp2arfm4p()
  initTestSuite;
end
function test_dtp2arfm4p_()
    x = load('sea.dat');                   % Load data 
   [dtp,u,tp] = dat2dtp([-2 2 32],x,0.2); % Discrete TP & rainflow filter 0.2 
   [ARFM,res] = dtp2arfm4p(dtp,32);      % Calculate asymmetric rainflow matrix 
   cmatplot(u,u,ARFM,3);                  % Plot rainflow matrix 
   colorbar;   
   % res  
 
   close all;
end
