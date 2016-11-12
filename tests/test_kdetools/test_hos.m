function test_suite=test_hos()
  initTestSuite;
end
function test_hos_()
   data = rndnorm(0, 1,20,1) 
  h = hos(data,'epan');
end
