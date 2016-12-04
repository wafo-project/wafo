function test_suite=test_varname()
  initTestSuite;
end
function test_varname_()
  x2 = 'test'; 
 assert(varname(x2), {'x2'}) 
 assert(varname(1), {''})
end
