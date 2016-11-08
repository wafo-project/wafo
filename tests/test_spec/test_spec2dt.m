function test_suite=test_spec2dt()
  initTestSuite;
end
function test_spec2dt_()
   S = jonswap; 
  assert(spec2dt(S), 1.04719755119660, 1e-10)
end
