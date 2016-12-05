function test_suite=test_getmfileinfo()
  initTestSuite;
end
function test_getmfileinfo_()
           % the getmfileinfo function 
 
  info1 = getmfileinfo('getmfileinfo.m', 'math','trivia'); 
  assert(info1.todo.comment{1}(1:30), ' Unable to extract/distinguish') 
  assert(info1.synopsis{1}, 'info1 = getmfileinfo(mfile,varargin)'); 
  assert(info1.subroutine.names{1}, 'removedoubleblanks'); 
  assert(info1.subroutine.synopsis{1}, 'str = removedoubleblanks(str)') 
  assert(info1.functions_called(2:5)', {'all','any','brackets','cell'});
end
