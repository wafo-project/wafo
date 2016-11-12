function test_suite=test_hste()
  initTestSuite;
end
function test_hste_()
    x  = rndnorm(0,1,50,1); 
   hs = hste(x,'gauss');
end
