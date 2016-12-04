function test_suite=test_ddewhite()
  initTestSuite;
end
function test_ddewhite_()
  
  s  = '    Testing testing 1-2-3   '; 
  cs = {s, ' deo, deo, \\n  '}; 
  assert(ddewhite(s), 'Testing testing 1-2-3') 
  assert(ddewhite(cs), {'Testing testing 1-2-3', 'deo, deo, \\n'})
end
