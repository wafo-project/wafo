function test_suite=test_hmns()
  initTestSuite;
end
function test_hmns_()
     data = rndnorm(0, 1,20,2) 
    h = hmns(data,'epan');
end
