function test_suite=test_dat2dtp()
  initTestSuite;
end
function test_dat2dtp_()
    x = load('sea.dat'); x=x(1:200,:);  
   [dtp,u,tp] = dat2dtp([-2 2 32],x,0.5); % Discrete TP & rainflow filter 0.5 
   plot(x(:,1), x(:,2), tp(:,1), tp(:,2)); 
 
   close all;
end
