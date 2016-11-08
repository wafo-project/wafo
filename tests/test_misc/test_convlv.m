function test_suite=test_convlv()
  initTestSuite;
end
function test_convlv_()
  dx = 2*pi/100;  
 x = linspace(0,2*pi-dx,100)'; 
 y = cos(x); 
 c = savgol(2,2,4,1);          % differentiation filter 
 yf = convlv(y,c)/dx;          % derivative 
 yf2 = convlv(y,c,[],'p')/dx;   
 semilogy(x,abs(sin(x)+yf),x,abs(sin(x)+yf2),'g')
end
