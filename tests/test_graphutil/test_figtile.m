function test_suite=test_figtile()
  initTestSuite;
end
function test_figtile_()
  for ix=1:10,figure(ix);end 
 figtile             % tile all open figures 
 figtile 1:3 5 7     % tile figure 1,2,3,5 and 7 
 figtile pairs2 1:10 % tile figure 1 to 10 two at a time 
 figtile pairs3 1:10 % tile figure 1 to 10 three at a time      
 
 close all;
end
