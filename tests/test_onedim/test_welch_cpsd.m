function test_suite=test_welch_cpsd()
  initTestSuite;
end
function test_welch_cpsd_()
   Fs = 400; x = 0:1/Fs:1; w1 = 2*pi*30; w2 = 2*w1; 
  y1 = cos(w1*x) + randn(size(x)); 
  y2 = cos(w2*x) + randn(size(x)); 
  [Pxy,fi] = welch_cpsd(y1,y2); 
  plot(fi,Pxy); 
 
  close all;
end
