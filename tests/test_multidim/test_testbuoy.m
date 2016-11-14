function test_suite=test_testbuoy()
  initTestSuite;
end
function test_testbuoy_()
    N=5000;dt=0.4;f0=0.1;th0=0;h=50;xypos=[0 0 0 1 1;0 0 0 4 0; 0 0 0 5 0]; 
   D = testbuoy(N,dt,3,f0,th0,0,0,h); 
   S = dat2dspec(D,xypos,h); 
   plotspec(S); 
 
   close all;
end
