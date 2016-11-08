function test_suite=test_jacobianest()
  initTestSuite;
end
function test_jacobianest_()
   xdata = (0:.1:1)'; 
  ydata = 1+2*exp(0.75*xdata); 
  fun = @(c) ((c(1)+c(2)*exp(c(3)*xdata)) - ydata).^2; 
 
  [jac,err] = jacobianest(fun,[1 2 0.75]);  
  assert(jac, zeros(size(jac)), 1e-12); 
  assert(err<1e-12);
end
