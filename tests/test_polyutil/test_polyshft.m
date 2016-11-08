function test_suite=test_polyshft()
  initTestSuite;
end
function test_polyshft_()
   py = [1 0]; 
  px = polyshft(py,0,5); 
  assert(polyval(px,[0 2.5 5]), [-1 0 1 ], eps);  % This is the same as the line below 
  assert(polyval(py,[-1 0 1 ]), [-1 0 1 ], eps);
end
