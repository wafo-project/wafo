function test_suite=test_chi2gof2()
  initTestSuite;
end
function test_chi2gof2_()
    xs  = rndray(1,500,1); 
   xs2 = rndnorm(0,1,100000,1);  
   p   = chi2gof2(pdfnorm(xs),pdfnorm(xs2));
end
