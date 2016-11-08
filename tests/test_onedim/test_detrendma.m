function test_suite=test_detrendma()
  initTestSuite;
end
function test_detrendma_()
   x = linspace(0,1,200)'; 
  y = exp(x)+cos(5*2*pi*x)+1e-1*randn(size(x)); 
  [y0, tr] = detrendma(y,20); 
  plot(x,y,x,y0,'r',x,exp(x),'k',x,tr,'m'); 
   
  close all;
end
