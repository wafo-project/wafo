function test_suite=test_comnsize()
  initTestSuite;
end
function test_comnsize_()
    A = rand(4,5);B = 2;C = rand(4,5); 
   [csiz,A1,B1,C1] = comnsize(A,B,C); 
   csiz = comnsize(A,1:2);
end
