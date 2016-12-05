function test_suite=test_freadtxt()
  initTestSuite;
end
function test_freadtxt_()
  txt = freadtxt(fullfile(waforoot, 'fileutil', 'freadtxt.m')); 
 assert(txt(1:27), 'function t = freadtxt(file)');
end
