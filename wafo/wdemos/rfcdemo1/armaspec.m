function [R]=armaspec(C,A,s2,n)

% ARMASPEC   Computes the spectral density for an AR- or ARMA-model.
%   The process is governed by the system equation
%     A(q) * x(t) = C(q) * sqrt(s2) * e(t) 
%
% S = armaspec(C,A,s2,n)
%
%
% S   = Spectral density. [f1 S1; f2 S2; ... fn Sn]
%       (Frequencies in row 1 and spectral density in row 2.)
%
% C   = Coefficients in C-polynomials. [1 c_1 ... c_nc]
% A   = Coefficients in A-polynomials. [1 a_1 ... a_na]
% s2  = Innovation variance.
% n   = Number of calculated values.
%
% Example: AR(2)-process.
%   S = armaspec(1,[1 1 0.9],1,500);
%   plot(S(:,1),S(:,2))
% Example: ARMA(4,2)-process.
%   S = armaspec([1 0.05 -0.88],[1 -2.06 1.64 -0.98 0.41],4.84e-6,500);
%   plot(S(:,1),S(:,2))

% Copyright (c) 1997 by Pär Johannesson
% Toolbox: Rainflow Cycles for Switching Processes V.1.0, 2-Oct-1997

[H,w]=freqz(C,A,n);
R=real(s2*H.*conj(H));
f=w/(2*pi);

R = [f R];
