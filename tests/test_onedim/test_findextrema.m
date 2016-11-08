function test_suite=test_findextrema()
  initTestSuite;
end
function test_findextrema_()
   t = linspace(0,7*pi,250); x = sin(t); 
  ind = findextrema(x); 
  plot(t,x,'.',t(ind),x(ind),'r.'); 
 
  assert(ind', [19,54,90,126,161,197,232], eps); 
  close all;
end
