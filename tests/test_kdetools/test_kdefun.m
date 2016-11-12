function test_suite=test_kdefun()
  initTestSuite;
end
function test_kdefun_()
    data = rndray(1,500,1); 
   x = linspace(sqrt(eps),5,55); 
   plotnorm((data).^(.5)); % gives a straight line => L2 = 0.5 reasonable 
   f = kdefun(data,{'L2',.5},x); 
   plot(x,f,x,pdfray(x,1),'r'); 
 
   close all;
end
