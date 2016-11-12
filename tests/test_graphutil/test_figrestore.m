function test_suite=test_figrestore()
  initTestSuite;
end
function test_figrestore_()
   for ix = 1:5, figure(ix); end 
  figrestore('all');   %Restores all unhidden figures 
  figrestore;          %same as figrestore('all') 
  figrestore(gcf);     %Restores the current figure 
  figrestore(3);       %Restores figure 3 
  figrestore([2 4]);   %Restores figures 2 and 4 
  % or alternativel 
  figrestore 2 4; 
 
  close all;
end
