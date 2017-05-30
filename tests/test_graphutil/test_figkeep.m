function test_suite=test_figkeep()
  initTestSuite;
end
function test_figkeep_()
   for ix = 1:10,figure(ix),end 
  figkeep 1:3  5 7  
  % or  
  figkeep([1:3  5 7])   
 
  close all;
end
