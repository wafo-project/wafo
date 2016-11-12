function test_suite=test_figpile()
  initTestSuite;
end
function test_figpile_()
  for ix = 1:10, figure(ix),end 
 figpile             % pile all open figures 
 figpile 1:3 5 7     % pile figure 1,2,3,5 and 7 
 
 close all;
end
