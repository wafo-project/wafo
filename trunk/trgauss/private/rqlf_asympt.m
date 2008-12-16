function [Cb0,Cb1,a0,a1]=rqlf_asympt(u,g,b,S22,S21)
% RQLF_ASYMPT Gives first two terms in an asymptotic expansion of the
% crossing rate of a quadratic+linear
%
%   CALL:  [Cb0,Cb1,a0,a1] = rqlf_asympt(u,g,b,S22,S12);
%
% A function that gives the two first terms
% in an asymptotic expansion of the rate of crossings of a
% quadratic+linear form according to Breitung's method


g=g(:);
b=b(:);

if length(g)~=length(b)
   error('The lengths of b and g should be the same.')
end
%if u<=0;
%   error('u should be a positive number.')
%end
n=size(S22,1);
xtp=mindist(g,b,u,20);
xJ=xtp(:);
b0=sqrt(sum(xJ.^2));
yJ=xJ/b0;
PJ=eye(n)-yJ*yJ';
dg=-b0*(b+2*(g.*xJ));
GJ=-2*b0^2*diag(g)/sqrt(sum(dg.^2));
V=S22+S21*S21;
sJ2=yJ'*V*yJ;
sJ=sqrt(sJ2);
a0=sJ*sqrt(1-yJ'*S21*(eye(n)+GJ)*S21*yJ/sJ2)/sqrt(det(eye(n)+ ...
						  PJ*GJ*PJ));
aJ=-1/sJ*PJ*(eye(n)+GJ)*S21*yJ;
WJ1=inv(eye(n)+PJ*GJ*PJ+aJ*aJ');
WJ2=inv(eye(n)+PJ*GJ*PJ);
vJ1=PJ*S21*PJ*GJ*yJ;
vJ2=PJ*GJ*yJ;
vJ3=PJ*GJ*S21*yJ;
vJ4=PJ*GJ*PJ*V*yJ;
KJ1=yJ'*S21*GJ*yJ;
KJ2=yJ'*GJ*yJ;
QJ1=PJ*GJ*PJ;
QJ2=PJ*S21*PJ*GJ*PJ;
QJ3=QJ2+.5*KJ1*QJ1-vJ2*vJ3';
QJ4=-(yJ'*V*PJ*GJ*yJ)*PJ*GJ*PJ-2*PJ*GJ*PJ*V*yJ*yJ'*GJ*PJ+PJ*GJ* ...
       PJ*V*PJ*GJ*PJ;
a11=1/2*holmquist1(WJ1,QJ4)/sJ2-1/2*holmquist1(WJ1,vJ4*vJ4')/sJ2^2+ ...
    1/2*holmquist2(WJ1,QJ3,QJ3)/sJ2+holmquist2(WJ1,aJ*vJ4',QJ3)/ ...
       sJ^3+1/2*holmquist2(WJ1,aJ*aJ',vJ4*vJ4')/sJ2^2+1/8* ...
holmquist3(WJ1,vJ2*vJ2',QJ1,QJ1)-1/2*holmquist2(WJ1,vJ2*vJ2',QJ1)-1/8*(KJ2+1)*holmquist2(WJ1,QJ1,QJ1)+...
1/2*holmquist2(WJ1,vJ2*vJ4',QJ1)/sJ2;
a12=(1/2*holmberg1(WJ2,aJ,vJ1,QJ1)+holmberg1(WJ2,aJ,vJ2,QJ2)-1/2*holmberg1(WJ2,aJ,vJ3,QJ1)+...
    KJ1*holmberg1(WJ2,aJ,vJ2,QJ1)-holmberg1(WJ2,aJ,vJ3,vJ2*vJ2')-1/2*KJ2*holmberg1(WJ2,aJ,vJ3,QJ1))/sJ-...
    1/2*holmberg2(WJ2,aJ,vJ2,QJ3,QJ1)/sJ-1/2*holmberg2(WJ2,aJ,aJ,vJ2*vJ2',QJ1)-...
    1/8*(KJ2+1)*holmberg2(WJ2,aJ,aJ,QJ1,QJ1)+1/8*holmberg3(WJ2,aJ,aJ,vJ2*vJ2',QJ1,QJ1);
[det(WJ1) det(WJ2) sJ a11 a12];
a1=sJ*(sqrt(det(WJ1)*a11)+sqrt(det(WJ2))*a12);
Cb0=a0*exp(-b0^2/2)/pi;
Cb1=a1*exp(-b0^2/2)/pi/b0^2;







