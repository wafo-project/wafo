function test_suite=test_hboot()
  initTestSuite;
end
function test_hboot_()
    data = rndnorm(0, 1,20,1); 
   [hs hvec score] = hboot(data,'epan'); 
   plot(hvec,score); 
 
   close all;
end
