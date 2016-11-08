function test_suite=test_derivest()
  initTestSuite;
end
function test_derivest_()
   opts = derivest('defaults'); % default options 
  opts.DerivativeOrder = 2; 
 
  %First derivative of exp(x), at x == 1 
   [d,e]  =derivest(@(x) exp(x),1); 
   [d2,e2]=derivest(@(x) exp(x),1,opts); 
   assert(d, exp(1), 1e-12); 
   assert(d2, exp(1), 1e-10); 
 
  %Third derivative of x.^3+x.^4, at x = [0,1] 
   d3 = derivest(@(x) x.^3 + x.^4,[0 1],'deriv',3); 
   assert(d3, [6,30], 1e-10);
end
