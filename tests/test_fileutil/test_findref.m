function test_suite=test_findref()
  initTestSuite;
end
function test_findref_()
    % (including all "math" functions) used in the 
   % findref function 
  ref = findref('findref',3:5,'math'); 
  assert(ref, strvcat('error', 'exist', 'floor'));
end
