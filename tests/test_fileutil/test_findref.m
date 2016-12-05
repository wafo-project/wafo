function test_suite=test_findref()
  initTestSuite;
end
function test_findref_()
    % (including all "math" functions) used in the 
   % findref function 
  [folder, root] = fileparts(waforoot);
  if strcmpi(root,'wafo'),
     ref = findref('findref',3:5,'math'); 
     assert(ref, strvcat('error', 'exist', 'floor'));
  end
end
