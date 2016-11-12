function test_suite=test_fcolorbar()
  initTestSuite;
end
function test_fcolorbar_()
  
 [x,y,z]=peaks; v=-10:2:8;  
 figure(1);  [cs,h]=contourf(x,y,z,v); clabel(cs,h); colorbar 
 figure(2);  [cs,h]=contourf(x,y,z,v); clabel(cs,h); fcolorbar(v); 
           % And not using contourspecification: 
 figure(3);  [cs,h]=contourf(x,y,z); clabel(cs,h); fcolorbar(cs); 
 
           % Not equidistant contours. 
 figure(4);  v=[-8 -4 -2 -1 0 1 2 4 8];  
 [cs,h]=contourf(x,y,z,v); clabel(cs,h); fcolorbar(v); 
 
 close all;
end
