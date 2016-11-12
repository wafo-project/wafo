function test_suite=test_kernelstats()
  initTestSuite;
end
function test_kernelstats_()
   [mu2,R]=kernelstats('triweight'); 
 
  assert(mu2, 0.111111111111111, eps) 
  assert(R, 0.815850815850816, eps)
end
