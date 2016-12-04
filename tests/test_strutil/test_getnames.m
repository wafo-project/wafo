function test_suite=test_getnames()
  initTestSuite;
end
function test_getnames_()
   names = getnames('yy = foo(x1.^2,a_1*c,flag)'); 
  assert(names, strvcat('yy', 'foo', 'x1', 'a_1', 'c', 'flag'))
end
