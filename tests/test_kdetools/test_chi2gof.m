function test_suite=test_chi2gof()
  initTestSuite;
end
function test_chi2gof_()
    xs = rndray(1,500,1); 
   x  = linspace(-7,7,201);  
   p  = chi2gof(pdfnorm(xs),pdfnorm(x));
end
