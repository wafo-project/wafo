function test_suite=test_dat2crossind()
  initTestSuite;
end
function test_dat2crossind_()
    t = linspace(0,7*pi,250);  
   x = sin(t); 
   [ind, Nc] = dat2crossind(x,0.75,'u'); 
   plot(t,x,'.',t(ind),x(ind),'o'); 
 
   assert(ind', [ 10,81,152,224], eps); 
   close all;
end
