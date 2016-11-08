function test_suite=test_iscomnsize()
  initTestSuite;
end
function test_iscomnsize_()
    A = rand(4,5);B = 2;C = rand(4,5); 
   [iscsize,A1,B1,C1] = iscomnsize(A,B,C); 
   iscsize = iscomnsize(A,1:2);
end
