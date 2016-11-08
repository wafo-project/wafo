function test_suite=test_plotlc()
  initTestSuite;
end
function test_plotlc_()
   x = load('sea.dat'); 
  lc = dat2lc(x,0.2);  
  plotlc(lc); 
 
  assert(length(lc), 269); 
  close all;
end
