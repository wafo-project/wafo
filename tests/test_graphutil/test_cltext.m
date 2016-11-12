function test_suite=test_cltext()
  initTestSuite;
end
function test_cltext_()
   z  = peaks; 
  cl = max(z(:))-range(z(:))*(.1:.1:.9); 
  contour(z,cl); cltext(cl); 
 
  data = rndray(1,2000,2);  
  f = kdebin(data,{'kernel','epan','L2',.5,'inc',128}); 
  contour(f.x{:},f.f,f.cl); cltext(f.pl,1); 
 
  close all;
end
