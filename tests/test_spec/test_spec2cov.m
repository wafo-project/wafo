function test_suite=test_spec2cov()
  initTestSuite;
end
function test_spec2cov_()
     S   = jonswap;  
    S.S([1:40 100:end]) = 0;    
    Nt  = length(S.S)-1;   
    R   = spec2cov(S,0,Nt); 
    win = parzen(2*Nt+1); 
    R.R = R.R.*win(Nt+1:end); 
    S1  = cov2spec(R); 
    R2  = spec2cov(S1); 
    figure(1) 
    plotspec(S),hold on, plotspec(S1,'r') 
    figure(2) 
    covplot(R), hold on, covplot(R2,[],[],'r') 
    figure(3) 
    semilogy(abs(R2.R-R.R)), hold on, 
    semilogy(abs(S1.S-S.S)+1e-7,'r')   
   
    close all
end
