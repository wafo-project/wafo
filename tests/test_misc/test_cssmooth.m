function test_suite=test_cssmooth()
  initTestSuite;
end
function test_cssmooth_()
   x = linspace(0,1).'; 
  y = exp(x)+1e-1*randn(size(x)); 
  pp = cssmooth(x,y,.9);  
  plot(x,y,x,cssmooth(x,y,.99,x,0,0.01),'g',x,ppval(pp,x),'k',x,exp(x),'r'); 
 
  close all;
end
