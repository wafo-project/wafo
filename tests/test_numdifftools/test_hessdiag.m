function test_suite=test_hessdiag()
  initTestSuite;
end
function test_hessdiag_()
   [HD,err] = hessdiag(@(x) x(1) + x(2)^2 + x(3)^3, [1 2 3]); 
  assert(HD, [ 0,2,18], 1e-12); 
  assert(err < 1e-12);
end
