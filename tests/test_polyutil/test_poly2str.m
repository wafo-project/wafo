function test_suite=test_poly2str()
  initTestSuite;
end
function test_poly2str_()
   assert(poly2str( [1 1 2], 's' ), 's^2 + s + 2');
end
