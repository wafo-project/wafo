function test_suite=test_hstt()
  initTestSuite;
end
function test_hstt_()
    x  = rndnorm(0,1,50,1); 
   hs = hstt(x,'gauss');
end
