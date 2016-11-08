function test_suite=test_dat2steep()
  initTestSuite;
end
function test_dat2steep_()
   dt = 0.4; 
  xs = spec2sdat(specinterp(jonswap,dt),6000);  
  rate=8; method=1; 
  [S,H] = dat2steep(xs,rate,method); 
  plot(S,H,'.'); 
  xlabel('Vcf [m/s]'); 
  ylabel('Hd [m]'); 
 
  close all;
end
