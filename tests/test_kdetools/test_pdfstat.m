function test_suite=test_pdfstat()
  initTestSuite;
end
function test_pdfstat_()
    x1=linspace(0,6)'; 
   par = {2 2  3 2.5 .8}; 
   [m, v] = momweib2d(par{:}); %  mean and covariance 
   [m, v] = momweib2d(par{:},'condon',2,'cvar',x1); 
   plot(x1,m,'r--', x1,sqrt(v),'g-'); % conditional mean and standard deviation. 
  
   [z1,z2] = rndweib2d(par{:},500,1); 
   f = kdebin([z1,z2]); 
   [M1,V1]=pdfstat(f,2,x1); 
 
   plot(x1,M1,'r--',x1,sqrt(V1),'k-'); 
   title(' Conditional mean and standard deviation'); 
   legend('E(x1|x2)','std(x1|x2)'); 
   xlabel('x2'); 
 
   close all; 
  % Does not work correctly yet!
end
