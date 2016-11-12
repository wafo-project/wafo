function test_suite=test_tallibing()
  initTestSuite;
end
function test_tallibing_()
   [x,y,z] = peaks(20);  
  epcolor(x,y,z);  
  tallibing(x,y,z); 
   % pcolor(x,y,z); shading interp;  
 
  close all;
end
