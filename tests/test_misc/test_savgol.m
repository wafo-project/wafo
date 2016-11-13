function test_suite=test_savgol()
  initTestSuite;
end
function test_savgol_()
    x = linspace(0,1); 
   y = exp(x)+1e-1*randn(size(x)); 
   nl=2; nr=2; m=2; ld=0; 
   c = savgol(nl,nr,m,ld,length(x)); 
   yf = real(ifft(fft(y).*fft(c))); % convolution with pronounced end effects 
   yf2 = convlv(y,c);               % convolution with less end effects 
   plot(x,y,x,yf,x,yf2); 
   ld =1; m =4; 
   c = savgol(nl,nr,m,ld);          % differentiation filter 
   dyf = convlv(y,c)*gamma(ld+1)/(x(2)-x(1))^ld; % Derivative 
   ix = nl+1:length(x)-nr;                       % for all these ix 
   semilogy(x(ix),abs(y(ix)-dyf(ix))); % Error of the derivative 
 
   close all;
end
