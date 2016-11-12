function test_suite=test_kdebin()
  initTestSuite;
end
function test_kdebin_()
   data = rndray(1,500,1); 
                               %Box-Cox transform data before estimation 
  f = kdebin(data,{'L2',.5,'inc',64});  
  pdfplot(f); 
                               %Non-parametric transformation 
  g   = cdf2tr(edf(data,'wdata',false),mean(data),std(data)); 
  opt = kdeoptset('L2',{g},'inc',64);   
  f1  = kdebin(data,opt); 
  hold on; pdfplot(f1,'r'); hold off; 
 
  close all;
end
