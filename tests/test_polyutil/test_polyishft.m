function test_suite=test_polyishft()
  initTestSuite;
end
function test_polyishft_()
   px = [1 0]; 
  py = polyishft(px,-5,5); 
  assert(polyval(px,[-5 0 5]), [-5 0 5], 1e-10);  % This is the same as the line below 
  assert(polyval(py,[-1 0 1 ]), [-5, 0, 5], 1e-10);
end
