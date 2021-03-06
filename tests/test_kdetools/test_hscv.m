function test_suite=test_hscv()
  initTestSuite;
end
function test_hscv_()
     % x = rndnorm(0, 1,20,1) 
  x = [-0.0233845632050972   0.9070186193622006;... 
        0.6529594866766634   1.3689145060433903;... 
        0.4477857310723146  -0.6311953712037597;... 
       -1.9256785038579962   0.5886257667993168;... 
       -0.5290011931824666  -0.3602090880229930]; 
  [hs hvec score] = hscv(x,'epan'); 
  plot(hvec,score); 
  assert(hscv(x,'epan'), [1.65372274837226, 1.54783142555932], 1e-10); 
  assert(hscv(x,'biwe'), [1.95910655435708, 1.83366086838792], 1e-10); 
  assert(hscv(x,'tria'), [1.81671456092952, 1.70038653180873], 1e-10); 
  assert(hscv(x,'triw'), [2.22465977863959, 2.08221016488133], 1e-10); 
  assert(hscv(x,'gaus'), [0.747004920174087, 0.699172634367494], 1e-10); 
  assert(hscv(x,'rect'), [1.29983118294192, 1.21660027647305], 1e-10); 
  assert(hscv(x,'lapl'), [0.552612179133844, 0.517227266693945], 1e-10); 
  assert(hscv(x,'logi'), [0.417582992650980, 0.390844281147180], 1e-10); 
 
  close all;
end
