function test_suite=test_ismatlab()
  initTestSuite;
end
function test_ismatlab_()
   assert(ismatlab, exist('matlabroot','builtin')~=0)
end
