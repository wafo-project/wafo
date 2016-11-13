function test_suite=test_hypgf()
  initTestSuite;
end
function test_hypgf_()
  x = linspace(-.99,.99)'; 
 [Sn1,err1] = hypgf(1,1,1,x); 
 plot(x,abs(Sn1-1./(1-x))+eps,'b',x,err1+eps,'r'),set(gca,'yscale','log'); 
 [Sn2,err2] = hypgf(.5,.5,1.5,x.^2); 
 plot(x,abs(x.*Sn2-asin(x))+eps,'b',... 
      x,abs(x.*err2)+eps,'r'); 
 set(gca,'yscale','log'); 
 
 close all;
end
