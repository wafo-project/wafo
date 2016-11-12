function test_suite=test_epcolor()
  initTestSuite;
end
function test_epcolor_()
  [x,y,z]=peaks(20);  
 epcolor(x,y,z); 
 
 close all;
end
