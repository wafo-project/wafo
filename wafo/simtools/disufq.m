function [svec, dvec] = disufq(A,w,kw,h,g,nmin,nmax)
%DISUFQ  Return difference- and sum-frequency effects.
%
%  CALL:  vec         = disufq(A,w,kw,h,g,nmin,nmax)
%         [svec,dvec] = disufq(A,w,kw,h,g,nmin,nmax)
%
% vec        = svec + dvec = 2nd order frequency component (size m X n).
% svec       = 2nd order sum frequency component           (size m X n).
% dvec       = 2nd order difference frequency component    (size m X n).
% A          = complex amplitudes (size m X n). 
% w          = vector with angular frequencies (w>=0)
% kw         = vector with wavenumbers (kw>=0)
% h          = water depth             (h >=0)
% g          = constant acceleration of gravity
% nmin       = minimum index where abs(A(:,nmin)) is 
%              greater than zero.
% nmax       = maximum index where abs(A(:,nmax)) is 
%              greater than zero.
%
% DISUFQ returns the summation of difference frequency and sum 
% frequency effects in the vector vec = svec + dvec
% The 2'nd order contribution to the non-linear wave is then calculated by
% a simple 1D Fourier transform, real(FFT(vec)).
%
% Examples:
%  % Estimate non-linear component  
% xn = load('sea.dat');
% dT = xn(2,1)-xn(1,1);
% n  = length(xn);
% h  = 10000;
% aMax = max(abs(xn(:,2))); 
% wMax = sqrt(2*gravity/aMax);
% A    = ifft(xn(:,2));
% w    = linspace(0,pi/dT,n/2);
% kw   = w2k(w);  
% nmax = min(max(find(w<=wMax)));
% nmin = 2;  
% vec  = disufq(A.',w,kw,h,gravity,nmin,nmax);  
% x2   = real(fft(vec.'));
% plot(xn(:,1),xn(:,2),'b',xn(:,1),xn(:,2)-x2,'r',xn(:,1),x2,'g')
% legend('nonlinear','approx linear','approx 2nd order comp')  
%
% % Simulate non-linear waves
% S    = jonswap(5);
% xs   = spec2sdat(S,n);
% A    = ifft(xs(:,2));
% w    = linspace(0,S.w(end),n/2);
% kw   = w2k(w);  
% nmax = min(max(find(w<=wMax)));
% nmin = 2;  
% vec  = disufq(A.',w,kw,h,gravity,nmin,nmax);  
% x2s  = real(fft(vec.'));
% plot(xs(:,1),xs(:,2),'b',xs(:,1),xs(:,2)+x2s,'r',xs(:,1),x2s,'g')
% legend('linear','non-linear','2nd order comp')      
% 
% close all;
% 
% See also  spec2nlsdat, spec2linspec  


% Is an internal function to spec2nlsdat
% History  
% This is a MEX-file for MATLAB.
% by Per Andreas Brodtkorb 15.08.2001
% revised pab 14.03.2002, 01.05.2002 22.07.2002

%error(nargchk(7,7,nargin))
narginchk(7,7)
disp('This function is only available as a mex-compiled function.')
error('Compile disufq.c by using mex -O disufq.c and try again.')
