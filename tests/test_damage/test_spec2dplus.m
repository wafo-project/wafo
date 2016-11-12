function test_suite=test_spec2dplus()
  initTestSuite;
end
function test_spec2dplus_()
   S = oscspec; bet = 3:0.2:5; 
  dplus = spec2dplus(S,bet); 
  assert(dplus, [ 2.32262632058483, 2.93360907105837,3.72269540324765,... 
                  4.74521181972958, 6.07455506745236, 7.80833828497044,... 
                 10.07670105706864, 13.05356835488749, 16.97193573771180,... 
                 22.14466211791423, 28.99281030590514], 1e-10);
end
