function test_suite=test_figcycle()
  initTestSuite;
end
function test_figcycle_()
   for ix=1:4,figure(ix),contourf(peaks(30)),end   
  % figcycle(1:3)          %Cycle trough figure 1 to 3 
  % figcycle 1:3           %Cycle trough figure 1 to 3  
  % figcycle 1:3 maximize  %Cycle trough figure 1 to 3 with figs maximized    
  % figcycle          %Cycle through all figures one at a time 
  % figtile pairs2    
  % figcycle pairs2  % Cycle through all figures two at a time  
 
  close all;
end
