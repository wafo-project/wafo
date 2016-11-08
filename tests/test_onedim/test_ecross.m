function test_suite=test_ecross()
  initTestSuite;
end
function test_ecross_()
   t = linspace(0,7*pi,250); x = sin(t); 
  ind = findcross(x,0.75);   
  t0 = ecross(t,x,ind,0.75); 
  plot(t,x,'.',t(ind),x(ind),'r.',... 
       t, ones(size(t))*.75, t0,ones(size(t0))*0.75,'g.'); 
 
  assert(ind', [10, 26, 81, 98, 152, 169, 224, 240], eps); 
  close all;
end
