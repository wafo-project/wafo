function test_suite=test_strflushleft()
  initTestSuite;
end
function test_strflushleft_()
   lines = {'Hello world', 'Hello space'}; 
  assert(strflushleft(lines, 7), {'       Hello world'; ... 
                                  '       Hello space'})
end
