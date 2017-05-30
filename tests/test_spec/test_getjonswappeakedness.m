function test_suite=test_getjonswappeakedness()
  initTestSuite;
end
function test_getjonswappeakedness_()
  assert(getjonswappeakedness(7,[5, 11, 25]), [7, 2.38529836797459   1], 1e-10) 
 Hm0 = linspace(1,20); 
 Tp = Hm0; 
 [T,H] = meshgrid(Tp,Hm0); 
 gam = getjonswappeakedness(H,T); 
 contourf(Tp,Hm0,gam,1:7),fcolorbar(1:7) 
 
 Hm0 = 1:10;   
 Tp  = linspace(2,16); 
 [T,H] = meshgrid(Tp,Hm0); 
 gam =  getjonswappeakedness(H,T);  
 plot(Tp,gam) 
 xlabel('Tp [s]')   
 ylabel('Peakedness parameter') 
 
 close all
end
