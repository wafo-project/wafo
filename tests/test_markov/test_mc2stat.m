function test_suite=test_mc2stat()
  initTestSuite;
end
function test_mc2stat_()
    F = magic(5); 
   P = mat2tmat(F); 
   ro = mc2stat(P); 
   assert(ro, ones(1,5)*0.2, eps)
end
