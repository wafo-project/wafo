function test_suite=test_vtkwrite()
  initTestSuite;
end
function test_vtkwrite_()
     
  [x,y,z] = meshgrid(-2:.2:2, -2:.25:2, -2:.16:2); 
  v = x .* exp(-x.^2 - y.^2 - z.^2); 
  vtkwrite(x,y,z,v,'test.vtk')
end
