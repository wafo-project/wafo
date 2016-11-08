function test_suite=test_wdensity()
  initTestSuite;
end
function test_wdensity_()
   S = linspace(20,40)'; T = linspace(4,20)'; 
  [S1 T1]=meshgrid(S,T); 
  sc = contour(S,T,wdensity(S1,T1)); clabel(sc) 
  xlabel('Salinity'),ylabel('Temperature')
end
