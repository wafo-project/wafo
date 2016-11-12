function test_suite=test_hlscv()
  initTestSuite;
end
function test_hlscv_()
   % x = rndnorm(0, 1,5,2); 
  x = [-0.0233845632050972   0.9070186193622006;... 
        0.6529594866766634   1.3689145060433903;... 
        0.4477857310723146  -0.6311953712037597;... 
       -1.9256785038579962   0.5886257667993168;... 
       -0.5290011931824666  -0.3602090880229930]; 
  [hs hvec score] = hlscv(x,'epan'); 
  plot(hvec,score); 
  assert(hlscv(x,'epan'), 1.39391122836191, 1e-10); 
  assert(hlscv(x,'biwe'), 1.65131708223985, 1e-10); 
  assert(hlscv(x,'triw'), 1.87515002002791, 1e-10); 
  assert(hlscv(x,'tria'), 1.53129587634976, 1e-10); 
  assert(hlscv(x,'gaus'), 0.629645172927051, 1e-10); 
  assert(hlscv(x,'rect'), 1.02091726336704, 1e-10); 
  assert(hlscv(x,'lapl'), 0.465792904029649, 1e-10); 
  assert(hlscv(x,'logi'), 0.351977756127559, 1e-10); 
 
  close all;
end
