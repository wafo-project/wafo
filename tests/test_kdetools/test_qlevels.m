function test_suite=test_qlevels()
  initTestSuite;
end
function test_qlevels_()
    x   = linspace(-8,8,2001); 
   qls = qlevels(pdfnorm(x),[10:20:90 95 99 99.9],x); 
 % compared with the exact values 
   ql  = pdfnorm(invnorm((100-[10:20:90 95 99 99.9])/200));   
 
   close all;
end
