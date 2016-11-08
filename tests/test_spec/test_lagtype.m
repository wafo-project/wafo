function test_suite=test_lagtype()
  initTestSuite;
end
function test_lagtype_()
   R = spec2cov(jonswap); 
  assert(lagtype(R), 't')
end
