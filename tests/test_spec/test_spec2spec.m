function test_suite=test_spec2spec()
  initTestSuite;
end
function test_spec2spec_()
     S    = demospec('dir'); 
    Snew = spec2spec(S,'enc',pi/6,10);
end
