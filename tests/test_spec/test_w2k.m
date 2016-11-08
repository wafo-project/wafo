function test_suite=test_w2k()
  initTestSuite;
end
function test_w2k_()
   w = linspace(0,3); 
  plot(w,w2k(w)) 
  assert(w2k([2.21430506782230, 3.13150025814577,... 
             4.42861013564459, 6.26300051629153]), [.5, 1, 2, 4], 1e-10) 
 
 close all
end
