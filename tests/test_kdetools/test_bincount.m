function test_suite=test_bincount()
  initTestSuite;
end
function test_bincount_()
   N  = 500; dx = 0.2; 
  f  = rndray(1,N,1); 
  ix = floor(f/dx)+1; 
  [len,bin] = bincount(ix);  
  plot((bin-.5)*dx,len/N/dx,'.'); % 1D probability density plot 
  bar((bin-.5)*dx,len/N/dx);      % 1D probability density plot 
  bar((bin-.5)*dx,len);           % 1D Histogram 
 
  close all;
end
