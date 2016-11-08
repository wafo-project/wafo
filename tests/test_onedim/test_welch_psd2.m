function test_suite=test_welch_psd2()
  initTestSuite;
end
function test_welch_psd2_()
  x = load('sea.dat'); 
 Fs = 1/diff(x(1:2,1)); 
 [Si,fi] = welch_psd2(x(:,2),[],Fs); 
 plot(fi, Si); 
 
 close all;
end
