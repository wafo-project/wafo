function [out]=ewwdir(omega,theta,omegat,thetat,h)
% EWWDIR Computes values of the quadratic transfer function E, for quadratic sea 
%
%  CALL:  [Eww]= ewwdir(w,th,wt,tht,h);
%  
%     Eww  = a matrix with the quadratic transfer function E(w,th,wt,tht). 
%     w,wt = two equally long vectors with angular frequencies.
%     w,wt = two equally long vectors with angular frequencies.
%     h    = water depth (default 5000 [m]).
%
%   Function uses w2k and is used in  dirsp2chitwo

%
%----------------------------------------------------------------------
% References: Marc Prevosto "Statistics of wave crests from second
% order irregular wave 3D models"
%
% Reduces to E(w,wt) from Eq.(6) in R. Butler, U. Machado, I. Rychlik (2002) 
% if th, tht are constant - longcrested sea eww.m.
% By I.R 22.10.04

g=gravity;
eps0=0.000001;
if ((length(omega)~=length(omegat))||(length(theta)~=length(thetat))||(length(thetat)~=length(omegat)))
   error('error in input to eww_new')    
end

if nargin<5 || h<=0
  h=5000;
end

[w wt]=meshgrid(omega,omegat);
[th tht]=meshgrid(theta,thetat);
wpl=w+wt;


ind=find(abs(w.*wt)<eps0);
ind1=find(abs(wpl)<eps0);
wpl(ind)=1;
w(ind)=1;
wt(ind)=1;
	 
   kw=w2k(w,[],h,g);
   kwx=kw.*cos(th);   kwy=kw.*sin(th);
   kwt=w2k(wt,[],h,g);
   kwtx=kwt.*cos(tht);   kwty=kwt.*sin(tht);
   kk=sqrt((kwx+kwtx).^2+(kwy+kwty).^2); kkh=g*kk.*tanh(kk*h);
   Dkwkwt=(2*wpl.*(g^2*((kwx.*kwtx)+(kwy.*kwty))-(w.*wt).^2)+g^2*((kw.^2.*wt)+(kwt.^2.*w))-wt.*w.*(wt.^3+w.^3))...
   ./(2*w.*wt.*(wpl.^2-kkh));
   Dkwkwt(ind1)=0.;
   out=(1/2/g)*(-g^2*((kwx.*kwtx)+(kwy.*kwty))./(w.*wt)+w.^2+wt.^2+w.*wt+2*wpl.*Dkwkwt);
   out(ind)=0;
   