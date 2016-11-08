function test_suite=test_findcross()
  initTestSuite;
end
function test_findcross_()
   v = 0.75; 
  t = linspace(0,7*pi,250); x = sin(t); 
  ind = findcross(x,v); 
  plot(t,x,'.',t(ind),x(ind),'r.', t, ones(size(t))*v); 
 
  assert(ind', [ 10,26,81,98,152,169, 224, 240], eps); 
  close all;
end
