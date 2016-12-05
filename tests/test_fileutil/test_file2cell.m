function test_suite=test_file2cell()
  initTestSuite;
end
function test_file2cell_()
  
 list = file2cell('file2cell.m'); 
 assert(list{1},'function list = file2cell(file, opt)');
end
