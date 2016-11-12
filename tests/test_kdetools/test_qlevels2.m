function test_suite=test_qlevels2()
  initTestSuite;
end
function test_qlevels2_()
    xs  = rndnorm(0,1,100000,1); 
   qls = qlevels2(pdfnorm(xs),[10:20:90 95 99 99.9]); 
            % compared with the exact values 
   ql  = pdfnorm(invnorm((100-[10:20:90 95 99 99.9])/200)); 
 
 % Finding the median of xs: 
   ql  = qlevels2(xs,50);
end
