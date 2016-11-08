function test_suite=test_spwaveplot()
  initTestSuite;
end
function test_spwaveplot_()
   x = load('sea.dat'); x1 = x(1:500,:); 
  spwaveplot(x1,[6:8 12:17]); 
 
  close all;
end
