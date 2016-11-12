function test_suite=test_plotfill()
  initTestSuite;
end
function test_plotfill_()
   x = linspace(0,10).'; 
  y = [sin(x)+1,sin(x)-1]; 
  plotfill(x,y,'g','facealpha',0.1) 
 
  close all;
end
