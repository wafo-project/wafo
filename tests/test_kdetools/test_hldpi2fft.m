function test_suite=test_hldpi2fft()
  initTestSuite;
end
function test_hldpi2fft_()
    x  = rndnorm(0,1,50,2); 
   hs = hldpi2fft(x,'gauss',1);
end
