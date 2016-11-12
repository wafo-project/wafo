function test_suite=test_gridcount()
  initTestSuite;
end
function test_gridcount_()
   N     = 500; 
  data  = rndray(1,N,1); 
  x = linspace(0,max(data)+1,50).';   
  dx = x(2)-x(1);   
  c = gridcount(data,x); 
  plot(x,c,'.')   % 1D histogram 
  plot(x,c/dx/N)  % 1D probability density plot 
  % trapz(x,c/dx/N)    
 
  close all;
end
