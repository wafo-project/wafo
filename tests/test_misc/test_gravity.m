function test_suite=test_gravity()
  initTestSuite;
end
function test_gravity_()
     assert(gravity, 9.80629386676700, 1e-10);
end
