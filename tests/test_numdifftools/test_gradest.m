function test_suite=test_gradest()
  initTestSuite;
end
function test_gradest_()
   [grad,err] = gradest(@(x) sum(x.^2),[1 2 3]);  
  assert(grad, [ 2,4, 6], 1e-12); 
  assert(err < 1e-12); 
 
  %At [x,y] = [1,1], compute the numerical gradient 
  %of the function sin(x-y) + y*exp(x) 
 
  z = @(xy) sin(-diff(xy)) + xy(2)*exp(xy(1)); 
 
  [grad2,err2 ] = gradest(z,[1 1]); 
  assert(grad2, [3.71828182845911,   1.71828182845906], 1e-12); 
  assert(err2 < 1e-12); 
 
  %At the global minimizer (1,1) of the Rosenbrock function, 
  %compute the gradient. It should be essentially zero. 
 
  rosen = @(x) (1-x(1)).^2 + 105*(x(2)-x(1).^2).^2; 
  [grad3,err3] = gradest(rosen,[1 1]); 
  assert(grad3, [0, 0], 1e-12); 
  assert(err3<1e-12);
end
