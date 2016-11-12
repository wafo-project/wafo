function test_suite=test_hns()
  initTestSuite;
end
function test_hns_()
     data = rndnorm(0, 1,20,1) 
    h = hns(data,'epan');
end
