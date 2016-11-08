function test_suite=test_loaddata()
  initTestSuite;
end
function test_loaddata_()
  
   x = loaddata('sea.dat');  
   assert(size(x), [9524, 2]) 
   assert(x(1:3,2)', [ -1.2004945, -1.0904945, -0.79049454], 1e-6) 
   assert(x(1:3,1)', [0.05, 0.3, 0.55], 1e-6)
end
