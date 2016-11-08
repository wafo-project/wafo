function test_suite=test_specplot()
  initTestSuite;
end
function test_specplot_()
   S = demospec('dir'); S2 = mkdspec(jonswap,spreading); 
  specplot(S,2), hold on 
  specplot(S,3,'g')  % Same as previous fig. due to frequency independent spreading 
  specplot(S2,2,'r') % Not the same as previous figs. due to frequency dependent spreading 
  specplot(S2,3,'m') 
  % transform from angular frequency and radians to frequency and degrees 
  Sf = ttspec(S,'f','d'); clf 
  specplot(Sf,2), 
 
  close all
end
