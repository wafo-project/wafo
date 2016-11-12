function test_suite=test_figstack()
  initTestSuite;
end
function test_figstack_()
  figstack             % stack all open figures 
 figstack 1:3 5 7     % stack figure 1,2,3,5 and 7 
 
 close all;
end
