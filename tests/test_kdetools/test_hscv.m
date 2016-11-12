function test_suite=test_hscv()
  initTestSuite;
end
function test_hscv_()
    data = rndnorm(0,1,20,1) 
   [hs hvec score] = hscv(data,'epan'); 
   plot(hvec,score)  
 
   close all;
end
