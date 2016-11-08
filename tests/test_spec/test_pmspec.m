function test_suite=test_pmspec()
  initTestSuite;
end
function test_pmspec_()
  S = pmspec(1.5,[6.5 10]); plotspec(S)
end
