function test_suite=test_clevels()
  initTestSuite;
end
function test_clevels_()
   
 c = contour(peaks); 
 [cl, c1] = clevels(c); 
 plot(c1(1,:),c1(2,:)); 
 cltext(cl) 
 
 close all;
end
