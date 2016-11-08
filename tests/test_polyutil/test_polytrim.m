function test_suite=test_polytrim()
  initTestSuite;
end
function test_polytrim_()
    assert(polytrim([0,0, 1,1,1]), [1,1,1], 1e-10);
end
