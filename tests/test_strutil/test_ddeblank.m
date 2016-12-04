function test_suite=test_ddeblank()
  initTestSuite;
end
function test_ddeblank_()
   s  = '    Testing testing 1-2-3   '; 
  cs = {s, ' deo, deo, \\n  '}; 
  assert(ddeblank(s), 'Testing testing 1-2-3') 
  assert(ddeblank(cs), {'Testing testing 1-2-3', 'deo, deo, \\n'})
end
