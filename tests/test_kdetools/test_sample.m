function test_suite=test_sample()
  initTestSuite;
end
function test_sample_()
      data = rndnorm(0,1,500,3); 
     s    = sample(data,100,0);
end
