function test_suite=test_gaussq()
  initTestSuite;
end
function test_gaussq_()
            % b) integration of x^2*exp(-x) from zero to infinity  
           % c) integrate humps            from 0 to 2 and from 1 to 4 
 
  A = [0 1]; B = [2,4]; 
  [val1,err1] = gaussq('(x.^2)',A,B);                 % a) 
  [val2,err2] = gaussq('(1)',0,inf,[1e-3 3],[],2);    % b) 
  [val2,err2] = gaussq('(x.^2)',0,inf,[1e-3 3],[],0); % b)   
  [val3,err3] = gaussq(@humps,A,B);                   % c) 
  assert(val1, [2.66666666666667, 21], 1e-10) 
  assert(val2, 2, 1e-10) 
  assert(val3, [34.92621394352225, 5.76237544601305], 1e-7)
end
