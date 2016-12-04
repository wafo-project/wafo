function test_suite=test_strcenter()
  initTestSuite;
end
function test_strcenter_()
  
  lines = {'Hello world', 'Hello space'}; 
  assert(strcenter(lines, 25), {'       Hello world'; ... 
                                '       Hello space'})
end
