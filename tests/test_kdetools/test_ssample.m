function test_suite=test_ssample()
  initTestSuite;
end
function test_ssample_()
      data = rndnorm(0,1,500,1); 
     s = ssample(data,100,'kernel','gauss'); 
     f = kdebin(s);   
     pdfplot(f); hold on; 
     x = linspace(-5,5); 
     plot(x,pdfnorm(x),'r'); 
 
     close all;
end
