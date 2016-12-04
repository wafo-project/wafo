function test_suite=test_strflushright()
  initTestSuite;
end
function test_strflushright_()
   lines = {'Hello world', 'Hello space'}; 
  assert(strflushright(lines, 18), {'       Hello world'; ... 
                                    '       Hello space'})
end
