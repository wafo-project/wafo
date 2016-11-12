function test_suite=test_figmaximize()
  initTestSuite;
end
function test_figmaximize_()
   for ix = 1:5, figure(ix),end 
  figmaximize('all')   %Maximizes all unhidden figures 
  figmaximize          %same as figmaximize('all') 
  figmaximize(gcf)     %Maximizes the current figure  
  figmaximize(3)       %Maximizes figure 3 
  figmaximize([2 4])   %Maximizes figures 2 and 4 
  figmaximize(gcf,'left')  %Windows taskbar at left of screen  
 %or alternatively 
  figmaximize 2 4   
 
  close all;
end
