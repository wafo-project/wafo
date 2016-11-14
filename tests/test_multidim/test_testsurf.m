function test_suite=test_testsurf()
  initTestSuite;
end
function test_testsurf_()
    N=5000;dt=0.4;f0=0.1;th0=0;h=50;xypos = [0 0 0 1 1;0 40 0 1 1; 20 20 0 1 1]; 
   D = testsurf(N,dt,3,f0,th0,xypos(:,1),xypos(:,2),h); 
   S = dat2dspec(D,xypos,h); 
   plotspec(S); 
 
   close all;
end
