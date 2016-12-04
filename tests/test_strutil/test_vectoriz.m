function test_suite=test_vectoriz()
  initTestSuite;
end
function test_vectoriz_()
  assert(vectoriz('x^2 + 1/x + x*y'), 'x.^2 + 1./x + x.*y')
end
