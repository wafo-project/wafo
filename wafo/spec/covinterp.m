function Rout=covinterp(R,t,vari)
% COVINTERP Interpolation of covariance function and derivatives
% 
% CALL: Rout = covinterp(Rin,t,vari);
%
%       Rout = covariance structure with new grid
%       Rin  = covariance structure
%       t    = vector of new grid points 
%              OR t=[dt N] (used if length(t)=2)
%              OR t=dt, with N such that not interpolation out of range
%                 (used if length(t)=1)
%              (default no change of input)
%       vari = 'x', 'y' OR 't', dimension to interpolate in 
%              (default 't', or if no t-variable then 'x')
%
% Interpolates all the matrices in the input covariance structure  
% (i.e. covariance function and all given derivatives)
% w.r.t dimension given by 'vari' such that the output is given at specified 
% points. If input t=[dt N], then output t=(0:N-1)*dt.
% If input t=dt, then N=length(Rin.t), and  output t=(0:N-1)*dt.
%      
% See also  spec2cov

% tested on: Matlab 7.0
% History: 
% revised pab feb 2007
% -replaced call to setfield and getfield with the equivalent syntax R.(fieldname), which is much faster.  
% revised by es 30.01.00   
% by es 13.10.1999  
  
if nargin < 2
  return
end
onedim=0;
Nin=length(R.R);
if numel(size(R.R))==Nin, % one-dim cvf
  onedim=1;
end
names=fieldnames(R);
if nargin<3
  if onedim
    ind=find(strcmp(names,'x')+strcmp(names,'t')); %options are 't' and 'x'
    if length(ind)>1
      if length(R.t)>1
	vari='t';
      else
	vari='x';
      end
    else
      vari=lower(names{ind});
    end
  else
    vari='t';
  end
end

tin=R.(vari);  %tin is either R.t or R.x
tin=tin(:)'; % make it row
if min(tin)>=0
  tin=[-tin(4:-1:2) tin];
end
if length(t)<=2
  if length(t)<2
    N=floor(tin(end)/t);
  else
    N=t(2);
  end
  t=(0:N-1)*t(1);
end
if t(end)>tin(end)
  disp('Warning: interpolation outside range, NaN in output')
end

ind=find(strncmp(names,'R',1)); 
d=ndims(R.R);
if d==3
  error('Three-dimensional interpolation not available yet')
end
Rout=R;
for j=1:length(ind)
  Y = R.(names{ind(j)});
  if onedim
    Y = Y(:)'; % make sure it is a row
    Y = [Y(:,4:-1:2) Y]';
  elseif strcmp(vari,'t')
    Y = [Y(:,4:-1:2) Y]';
  end
  if rem(length(names{ind(j)}),2)==1 %even order derivative
    Y = interp1(tin,Y,t,'cubic*')';
  else
    Y = interp1(tin,Y,t,'cubic*')';
  end
  if  strcmp(vari,'x') && isfield(R,'t')
    Y = Y';
  end
  Rout.(names{ind(j)}) = Y;
end
Rout.(vari) = t;
