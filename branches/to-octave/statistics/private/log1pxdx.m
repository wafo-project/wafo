function y = log1pxdx(x)
%LOG1PXDX Computes Log(1+x)/x
%

y = ones(size(x));
p = x+1==0;
k = (x~=0 & ~p);
y(k) = log1p(x(k))./x(k);
y(p) = realmax;
y(x==inf) = 0;