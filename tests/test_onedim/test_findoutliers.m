function test_suite=test_findoutliers()
  initTestSuite;
end
function test_findoutliers_()
  xx = load('sea.dat'); 
 dt = diff(xx(1:2,1)); 
 dcrit = 5*dt; 
 ddcrit = 9.81/2*dt*dt; 
 zcrit = 0; 
 [inds, indg] = findoutliers(xx,zcrit,dcrit,ddcrit); 
 waveplot(xx,'-',xx(inds,:),1,1,1); 
 
 assert(length(inds), 1152); 
 close all;
end
