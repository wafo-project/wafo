function test_suite=test_mayaplot()
  initTestSuite;
end
function test_mayaplot_()
   [x,y,z] = meshgrid(-2:.2:2, -2:.25:2, -2:.16:2); 
  v = x .* exp(-x.^2 - y.^2 - z.^2); 
%  mayaplot(x,y,z,v)
end
