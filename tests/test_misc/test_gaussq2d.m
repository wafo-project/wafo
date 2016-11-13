function test_suite=test_gaussq2d()
  initTestSuite;
end
function test_gaussq2d_()
  
   p1=2.; p2=0.5; 
   fun = '(p2*x.^2.*y+p1)'; 
   assert(gaussq2d(fun,[0 2],[1 4],[-1 1],[3 2],[],p1,p2),[10, 12], 1e-10);
end
