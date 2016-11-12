function test_suite=test_kde()
  initTestSuite;
end
function test_kde_()
    data = rndray(1,500,1); 
   x = linspace(sqrt(eps),5,55); 
   plotnorm((data).^(.5)) % gives a straight line => L2 = 0.5 reasonable 
   f = kde(data,{'L2',.5},x); 
   pdfplot(f) 
 
   close all;
end
