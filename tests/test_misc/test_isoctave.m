function test_suite=test_isoctave()
  initTestSuite;
end
function test_isoctave_()
   
  assert(isoctave, exist('octave_core_file_name', 'builtin')~=0)
end
