function test_suite=test_poly2hstr()
  initTestSuite;
end
function test_poly2hstr_()
   assert(poly2hstr( [1 1 2], 's' ), '(s + 1)*s + 2');
end
