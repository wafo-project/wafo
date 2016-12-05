function test_suite=test_plotspec()
  initTestSuite;
end
function test_plotspec_()
   S = demospec('dir');  
  S2 = mkdspec(jonswap, spreading); 
  plotspec(S,2); hold on; 
  plotspec(S,3,'g');  % Same as previous fig. due to frequency independent spreading 
  plotspec(S2,2,'r'); % Not the same as previous figs. due to frequency dependent spreading 
  plotspec(S2,3,'m'); 
  % transform from angular frequency and radians to frequency and degrees 
  Sf = ttspec(S,'f','d'); clf 
  plotspec(Sf,2); 
 
  close all;
end
