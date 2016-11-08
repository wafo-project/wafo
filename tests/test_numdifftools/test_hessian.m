function test_suite=test_hessian()
  initTestSuite;
end
function test_hessian_()
   %Rosenbrock function, minimized at [1,1] 
  rosen = @(x) (1-x(1)).^2 + 105*(x(2)-x(1).^2).^2; 
  [h,err] = hessian(rosen,[1 1]);   
  assert(h, [ 842 -420; -420, 210], 1e-10); 
  assert(err < 1e-9); 
 
  %cos(x-y), at (0,0) 
  opts = hessian('defaults'); 
  opts.RombergTerms = 4; 
  [h2, err2] = hessian(@(xy) cos(xy(1)-xy(2)),[0 0]); 
  assert(h2, [-1 1; 1 -1], 1e-6); 
  assert(err2 < 1e-5); 
 
  [h3, err3] = hessian(@(xy) cos(xy(1)-xy(2)),[0 0],opts); 
  assert(h3, [-1 1; 1 -1], 1e-6); 
  assert(err3 < 1e-6);
end
