function test_suite=test_mkspreading()
  initTestSuite;
end
function test_mkspreading_()
   options = mkspreading('defaults'); % get default options. 
  options.sp(1) = 10;                % Set spa = 10 
  D = mkspreading('cos2s',options); 
  w = linspace(0,3,257); 
  theta = linspace(-pi,pi,129); 
  contour(D(theta,w)) 
 
        % Make frequency dependent direction 
   options.theta0 = inline('pi/6*w'); 
   D2 = mkspreading('cos2s',options); 
  contour(D2(theta,w)) 
 
  % Plot all spreading functions 
 alltypes = {'cos2s','box','mises','poisson','sech2','wnormal'}; 
 for ix =1:length( alltypes) 
   D3 = mkspreading(alltypes{ix},options); 
   figure(ix) 
   contour(D3(theta,w)),title(alltypes{ix}) 
 end 
 
 close all
end
