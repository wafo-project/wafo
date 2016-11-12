function test_suite=test_black()
  initTestSuite;
end
function test_black_()
  t = linspace(0,4*2*pi)';  
 x = sin(t);  
 y = cos(t); 
 figure(1); plot(t,x,t,y);  
 figure(2); black; plot(t,x,t,y) 
 black off; % reset ColorOrder and LinstyleOrder to factory settings  
 
 close all;
end
