function test_suite=test_directionaldiff()
  initTestSuite;
end
function test_directionaldiff_()
   %At the global minimizer (1,1) of the Rosenbrock function, 
  %compute the directional derivative in the direction [1 2] 
 
  vec = [1,2]; 
  rosen = @(x) (1-x(1)).^2 + 105*(x(2)-x(1).^2).^2; 
  [dd,err] = directionaldiff(rosen,[1 1],vec); % d = 0 
  assert(dd, 0, 1e-10)
end
