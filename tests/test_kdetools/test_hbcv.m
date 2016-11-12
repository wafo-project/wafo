function test_suite=test_hbcv()
  initTestSuite;
end
function test_hbcv_()
   % data = rndnorm(0, 1,5,2); 
  data = [-0.0233845632050972   0.9070186193622006;... 
           0.6529594866766634   1.3689145060433903;... 
           0.4477857310723146  -0.6311953712037597;... 
          -1.9256785038579962   0.5886257667993168;... 
          -0.5290011931824666  -0.3602090880229930]; 
   [hs hvec score] = hbcv(data,'epan'); 
   plot(hvec,score); 
   assert(hs,  1.39391122836191, 1e-10) 
   assert(hbcv(data,'biwe'), 1.65131708223985, 1e-10) 
   assert(hbcv(data,'triw'), 1.87515002002791, 1e-10) 
   assert(hbcv(data,'tria'), 1.53129587634976, 1e-10) 
   assert(hbcv(data,'gaus'), 0.629645172927051, 1e-10) 
   assert(hbcv(data,'rect'), 1.09561852654024, 1e-10) 
   assert(hbcv(data,'lapl'),  0.137620630736033, 1e-10) 
   assert(hbcv(data,'logi'),  0.0879944390318899, 1e-10) 
 
   close all;
end
