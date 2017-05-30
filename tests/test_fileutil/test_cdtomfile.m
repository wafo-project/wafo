function test_suite=test_cdtomfile()
  initTestSuite;
end
function test_cdtomfile_()
   p0 = pwd(); 
  cdtomfile('cdtomfile'); 
  p1 = pwd(); 
  assert(p1, fullfile(waforoot, 'fileutil')); 
  cd(p0);
end
