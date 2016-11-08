function test_suite=test_dat2lc()
  initTestSuite;
end
function test_dat2lc_()
   x = load('sea.dat');  
  lc = dat2lc(x,0.2,1); 
  plotlc(lc); 
  plot(lc(:,1),lc(:,2)); 
 
  close all;
end
