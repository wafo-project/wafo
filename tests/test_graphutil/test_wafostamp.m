function test_suite=test_wafostamp()
  initTestSuite;
end
function test_wafostamp_()
    plot(sin(0:.1:3)); wafostamp('Sinus plot','(ER)'); hold on; 
   plot(sin((0:.1:3)+pi/2)); hold off; 
 
   close all;
end
