function test_suite=test_strim()
  initTestSuite;
end
function test_strim_()
   s  = '    Testing testing 1-2-3   '; 
  assert(strim(s), 'Testing testing 1-2-3')
end
