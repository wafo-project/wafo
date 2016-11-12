function test_suite=test_hldpi2()
  initTestSuite;
end
function test_hldpi2_()
    x  = rndnorm(0,1,50,2); 
   hs = hldpi2(x,'gauss',1);
end
