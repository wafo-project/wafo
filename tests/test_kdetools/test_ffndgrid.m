function test_suite=test_ffndgrid()
  initTestSuite;
end
function test_ffndgrid_()
  N = 500;D=2; sz = [N ,D ]; 
 x = randn(sz); z = ones(sz(1),1);  
 [nc, xv] = ffndgrid(x,z,'delta',-15,'method','sum');  % Histogram 
 pcolor(xv{:},nc)     % 
 [XV,YV]=meshgrid(xv{:}); 
 text(XV(:),YV(:),int2str(nc(:))) 
 dx = [diff(xv{1}(1:2)) diff(xv{2}(1:2))]; 
 contourf(xv{:}, nc/(N*prod(dx))); % 2-D probability density plot. 
 colorbar; 
 colormap jet; 
 
 close all;
end
