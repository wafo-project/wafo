function test_suite=test_hlscv()
  initTestSuite;
end
function test_hlscv_()
    data = rndnorm(0,1,20,1) 
   [hs hvec score] = hlscv(data,'epan'); 
   plot(hvec,score); 
 
   close all;
end
