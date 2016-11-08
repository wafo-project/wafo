function test_suite=test_lcplot()
  initTestSuite;
end
function test_lcplot_()
   x = load('sea.dat'); 
  lc = dat2lc(x,0.2);  
  lcplot(lc); 
 
  close all;
end
