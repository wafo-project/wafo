function test_suite=test_hldpi()
  initTestSuite;
end
function test_hldpi_()
    x  = rndnorm(0,1,50,1); 
   hs = hldpi(x,'gauss',1);
end
