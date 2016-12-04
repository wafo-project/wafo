function test_suite=test_str2cell()
  initTestSuite;
end
function test_str2cell_()
  
  s = ['Hello world' char(10) 'Hello space']; 
  assert(str2cell(s), {'Hello world'; 'Hello space'})
end
