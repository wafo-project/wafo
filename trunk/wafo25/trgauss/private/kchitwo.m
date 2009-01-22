function [K]=kchitwo(s,t,gam,bet,S12,S22)
% KCHITWO Computes the cumulant generating function formula for noncentral Chi2-process.
%
% CALL: [K]=kchitwo(s,t,gam,bet,S12,S22);
%
% Note that s,t must be constants (vectors are not allowed):
%


% Reference: "Five lectures on reliability applications of Rice's formula for 
%             level crossings" Rychlik (2003).
%
% tested on matlab 6.1
% by IR 12.11.03, changed name and added some checks and comments.
% rev. IR 24.11.03, changed name and added some checks and comments.

n=length(gam);
V=S22-S12'*S12;
Gam=diag(gam);
I=eye(n);
A=I-2*s*Gam-2*t*(Gam*S12'+S12*Gam)-4*t*t*Gam*V*Gam;
T=(s*I+t*S12+2*t*t*Gam*V)*bet';
K=-log(det(A))/2+t*t*bet*V*bet'/2+T'*inv(A)*T/2;


