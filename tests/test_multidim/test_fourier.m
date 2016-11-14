function test_suite=test_fourier()
  initTestSuite;
end
function test_fourier_()
   T = 2*pi;M=5; 
  t = linspace(0,4*T).'; x = sin(t); 
  [a,b] = fourier(t,x,T,M); 
 
  assert(a, zeros(5,1), 1e-10); 
  assert(b, [0,4,0,0,0]', 1e-10)
end
