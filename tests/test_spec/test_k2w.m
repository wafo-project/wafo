function test_suite=test_k2w()
  initTestSuite;
end
function test_k2w_()
  assert(k2w([.5, 1, 2, 4]), [2.21430506782230, 3.13150025814577,... 
                             4.42861013564459, 6.26300051629153], 1e-10)
end
