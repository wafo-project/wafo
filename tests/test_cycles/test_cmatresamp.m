function test_suite=test_cmatresamp()
  initTestSuite;
end
function test_cmatresamp_()
    F = [0   1   3   1;... 
        0   0   2   5;... 
        0   0   0   4;... 
        0   0   0   0]; 
   FF = cmatresamp(F); 
 
   assert(sum(FF(:)), sum(F(:)));
end
