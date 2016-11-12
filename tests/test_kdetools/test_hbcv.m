function test_suite=test_hbcv()
  initTestSuite;
end
function test_hbcv_()
    data = rndnorm(0, 1,20,1) 
   [hs hvec score] = hbcv(data,'epan'); 
   plot(hvec,score); 
 
   close all;
end
