function test_suite=test_vc2sign()
  initTestSuite;
end
function test_vc2sign_()
   x = load('sea.dat'); 
  [vcf,h] = dat2steep(x); % extract Crest front velocity and wave height.   
  hs = vc2sign(h); 
  plot(hs(:,1),hs(:,2)); 
  xlabel('Ratio'); ylabel('Ratio-significant wave height.'); 
  hs3 = vc2sign(h,1/3); 
  assert(hs3(2),  1.77039326808989, 1e-10);   % significant wave height   
  S   = dat2spec(x); 
  assert(spec2char(S,'Hm0'), 1.89191906170386, 1e-10); 
 
  close all;
end
