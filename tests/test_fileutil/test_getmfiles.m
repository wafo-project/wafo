function test_suite=test_getmfiles()
  initTestSuite;
end
function test_getmfiles_()
    
  mfiles = getmfiles(fullfile(waforoot, 'fileutil')); 
  names = {}; 
  for i=1:length(mfiles), 
   [folder, names{i}] = fileparts(mfiles{i}); 
  end 
  assert(names(1:3), {'Contents','bindiff','cdtomfile'});
end
