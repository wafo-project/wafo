function test_suite=test_parseoptions()
  initTestSuite;
end
function test_parseoptions_()
   defaultoptions = struct('test',[],'integrate',[] ); 
  options = parseoptions(defaultoptions,'int','yes'); 
  assert(options, struct('test',[], 'integrate', 'yes')) 
  opt = defaultoptions; 
  opt.test = 'yes';   
  options2 = parseoptions(defaultoptions,'int','yes',opt); 
  assert(options2, struct('test','yes', 'integrate', 'yes'))
end
