function test_suite=test_getsubdirs()
  initTestSuite;
end
function test_getsubdirs_()
   d = getsubdirs(waforoot, 3); 
  names = {}; 
  for i=1:length(d), 
   [folder, names{i}] = fileparts(d{i}); 
  end 
  assert(names(2:4), {'@data_1d', '@data_2d' '@data_3d'})
end
