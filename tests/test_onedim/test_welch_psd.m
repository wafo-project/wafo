function test_suite=test_welch_psd()
  initTestSuite;
end
function test_welch_psd_()
     [b,a] = cheby1(4,3,[0.2, 0.4]);    % define noise colour 
    y = filter(b,a,randn(2^12,1)); 
    welch_psd(y);                      % estimate noise colour 
    opt = welch_psd('defaults'); 
    opt.nfft = 128; 
    opt.Fs = 3; 
    welch_psd(y,opt); 
    welch_psd(y,'nfft',128,'Fs',3);   % alternatively 
 
    close all;
end
