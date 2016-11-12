function test_suite=test_deriv2()
  initTestSuite;
end
function test_deriv2_()
  
    k42 = deriv2(0,0,'42'); 
    assert(k42, -0.477464829275686)
end
