function test_suite=test_tallibing3()
  initTestSuite;
end
function test_tallibing3_()
   [x,y,z] = peaks(20);  
  surf(x,y,z); tallibing3(x,y,z,z); 
 
  close all;
end
