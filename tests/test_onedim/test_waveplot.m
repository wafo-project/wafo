function test_suite=test_waveplot()
  initTestSuite;
end
function test_waveplot_()
  % Plot x1 with red lines and mark troughs and crests with  
 % blue circles. 
   x = load('sea.dat'); 
   x1 = x(1:150,:); 
   waveplot(x1,'r-','bo'); 
 
   close all;
end
