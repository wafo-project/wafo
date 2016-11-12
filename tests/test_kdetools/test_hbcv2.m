function test_suite=test_hbcv2()
  initTestSuite;
end
function test_hbcv2_()
     data = rndnorm(0, 1,20,2); 
    hs = hbcv2(data,'gaus');
end
