function test_suite=test_createcov()
  initTestSuite;
end
function test_createcov_()
   R = createcov(1,'tx');  
  R.date = ''; 
  trueR = struct('R',[], 'x', [], 't', [], 'h', inf, 'tr', [], 'phi', 0, ... 
       'type', 'none', 'norm', [], 'Rx', [], 'Rt', [], 'note', [], 'date', ''); 
  assert(fieldnames(R), fieldnames(trueR)) 
  assert(struct2cell(R), struct2cell(trueR))
end
