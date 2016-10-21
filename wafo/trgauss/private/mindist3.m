function [x0,pmx]=mindist3(g,b,u)
% MINDIST3 Finds point of minimal distance to the origin on the surface b'*x+g'*x.^2=u.
%
%   CALL:   [x0,pmx] = mindist3(g,b,u);
%
% Returns a point of minimal distance to the origin on the surface
% b'*x+g'*x.^2=u

g=g(:);
b=b(:);
d0=length(g);
if length(b)~=d0
   error('g and b should have the same length')
end
% If u==0, x0=0 is the unique solution
if u==0
   x0=zeros(d0,1);
   pmx=zeros(d0,1);
   return
end

nulldim=(g==0)&(b==0);% If b(i)=g(i)=0, x0(i)=0

if sum(nulldim)==d0% If all b's and g's are zeros and u=0, only zeros are returned
   if u==0
      x0=zeros(d0,1);
      pmx=zeros(d0,1);
      return
   else
      error('u is out of range')
   end
end

g(nulldim)=[];
b(nulldim)=[];
d=length(b);
if (all(g>0) || all(g<0)) && u==-.25*sum(b.^2./g)
   x00=-.5*b./g;
   pmx=zeros(d,1);
   x0=zeros(d0,1);
   x0(~nulldim)=x00;
end
if all(g>0) && u<-.25*sum(b.^2./g)
   error('u out of range')
end
if all(g<0) && u>-.25*sum(b.^2./g)
   error('u out of range')
end
pmx0=zeros(d,1);
x00=[inf;zeros(d-1,1)];
I0=b==0;
I1=b~=0;
gI0=g(I0);
gI0unique=unique(gI0);
for k=1:length(gI0unique)
   if all(g(I1)-gI0unique(k))
      xp=zeros(d,1);
      xp(I1)=.5*b(I1)./(gI0unique(k)-g(I1));
      up=u-sum(b.*xp+g.*xp.^2,1);
      if sign(up)*sign(gI0unique(k))>=0
	 xp(min(find(g==gI0unique(k))))=sqrt(up/gI0unique(k));
	 if sum(xp.^2)<sum(x00.^2)
	    x00=xp;
	    pmx0=zeros(d,1);
	    pmx0(g==gI0unique(k))=1;
	 end 
      end
   end 
end
% Construction of the rational function:
bet=b(I1);
gam=g(I1);
if sum(I1)>0
   R=[bet(1)^2*[0 2 -gam(1)];1 -2*gam(1) gam(1)^2];
   for k=2:sum(I1);
      R=rfadd(R,[bet(k)^2*[0 2 -gam(k)];1 -2*gam(k) gam(k)^2]);
   end
   lam=roots(R(1,:)-4*u*R(2,:));
   lam=real(lam(abs(real(lam))>(1e6)*abs(imag(lam))));
   x1=zeros(d,1);
   for k=1:length(lam)
      x1(I1)=.5*bet./(lam(k)-gam);
      if sum(x1.^2)<sum(x00.^2);
	 x00=x1;
	 pmx0=zeros(size(pmx0));
      end
   end
end


x0=zeros(d0,1);
x0(~nulldim)=x00;
pmx=zeros(d0,1);
pmx(~nulldim)=pmx0;
