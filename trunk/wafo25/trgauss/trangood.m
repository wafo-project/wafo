function f = trangood(ff,nmin,mini,maxi,nmax)
%TRANGOOD Makes a transformation that is suitable for efficient transforms.
%
%  CALL:  f = trangood(ff,nmin,mini,maxi,nmax);
%
%        f    = the good transform function, [X f(X)].
%        ff   = the input transform function, [X f(X)]. 
%        nmin = the minimum number of points in the good transform.
%               (Default  size(ff,1))
%        mini = the minimum data value to transform. 
%               (Default  min(ff(:,1)))
%        maxi = the maximum data value to transform.
%               (Default  max(ff(:,1)))
%        nmax = then maximum number of points in the good transform
%               (default inf)
%
% TRANGOOD interpolates  ff  linearly and optionally
%  extrapolates  ff  linearly outside the range of  ff(:,1)
%  with X uniformly spaced.
%
% See also   tranproc, interp1q


% History:
% revised pab 07.02.2001
% - added nmax
% - replaced x = linspace(f(1,1),f(nf,1),nmin)' with x = (f(1,1):df:f(nf,1))'
% revised pab 07.01.2001
%  -fixed a bug: x = linspace(f(1,1),f(nf,1),nmin)' is safer than using
%    x=f(1,1)+(0:nmin-1)'/(nmin-1)*(f(nf,1)-f(1,1));  with interp1q
% revised pab 12.11.2000
%  - updated header info: A more detailed description of what TRANGOOD does.
%  - changed interpolation with a call to interp1q which is much faster
%  added nargchk and isempty(maxi),....isempty(nmin)
% by ???

error(nargchk(1,5,nargin))
if (size(ff,2)~=2)
  error('ff  must be a two column matrix.')
end
if (size(ff,1)<2)
  error('ff  must have at least two rows.')
end

[f,i] = sort(ff(:,1));
f     = [f ff(i,2)];
clear i;
df    = diff(f(:,1));
if ( any(df<=0)), %eps 
  error('Duplicate x-values in  ff  not allowed.')
end

nf = size(f,1);
if (nargin<5)||isempty(nmax),  nmax = inf; end
if (nargin<4)||isempty(maxi),  maxi = f(nf,1); end
if (nargin<3)||isempty(mini),  mini = f(1,1);  end
if (nargin<2)||isempty(nmin),  nmin = nf;      end
if (nmin<2),    nmin = 2;  end
if (nmax<2),    nmax = 2;  end

if ( (nf<nmin) ||(nmax<nf) || any(abs(diff(df))>10*eps*(f(nf,1)-f(1,1))) )
  % pab 07.01.2001: Always choose the stepsize df so that 
  % it is an exactly representable number.
  % This is important when calculating numerical derivatives and is 
  % accomplished by the following.
  ni = min(nmin,nmax);
  if ni==2,
      x = f([1,nf],1);
  else
    df = (f(nf,1)-f(1,1))/(ni-1);
    df = donothing(df+2)-2;
    x = (f(1,1):df:f(nf,1)).';
  end
  % New call pab 11.11.2000: This is much quicker
  f = [ x interp1q(f(:,1),f(:,2),x)];     
  %f = [ x interp1(f(:,1),f(:,2),x,'linear')];  
end
% f(:,1) is now uniformly spaced
df = f(2,1)-f(1,1);

% Extrapolate linearly outside the range of ff
%---------------------------------------------- 
if (mini<f(1,1)),
  f1 = df*(floor((mini-f(1,1))/df):1:-1)';
  f2 = f(1,2)+f1*(f(2,2)-f(1,2))/(f(2,1)-f(1,1));
  f  = [f1+f(1,1) f2;f];
end
n = size(f,1);
if (maxi>f(n,1))
  f1 = df*(1:1:ceil((maxi-f(n,1))/df))';
  f2 = f(n,2)+f1*(f(n,2)-f(n-1,2))/(f(n,1)-f(n-1,1));
  f  = [f;f1+f(n,1) f2];
end

return
function y=donothing(x)
  y=x;
return

% Old call: (Saved it just in case...)
%  x=f(1,1)+(0:nmin-1)'/(nmin-1)*(f(nf,1)-f(1,1));
%  % "y=interp1(f(:,1),f(:,2),x,'linear');" is slow:
%  % Use the fact that transforms often are "continuous" to 
%  % find y_k=f(x_k) incrementally.
%  y=zeros(nmin,1);
%  i = 1;
%  for k=1:nmin
%    xx=x(k);
%    if (xx>f(i+1,1))
%      while (xx>f(i+1,1))
%        i=i+1;
%        if (i>=nf), i=nf-1; break;  end;
%      end
%    else % xx<=f(i+1,1)
%      while (xx<=f(i,1))
%        i=i-1;
%        if (i<1),   i=1;    break;  end;
%      end
%    end
%    x0 = f(i,1); x1 = f(i+1,1);
%    y0 = f(i,2); y1 = f(i+1,2);
%    y(k) = (xx-x0)*(y1-y0)/(x1-x0)+y0;
%  end
%  f=[x y];
%  clear x y;
% 
% 
