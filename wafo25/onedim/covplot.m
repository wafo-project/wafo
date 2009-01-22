function covplot(rn,L,plotflag,lintype)
%COVPLOT Plots the auto covariance function (ACF) 1D or 2D. 
%
% CALL:  covplot(r,L,plotflag,lintype)
%
%         r =  covariance structure
%  for one-dim:
%         L = Maximum lag for which the ACF is plotted.
%     plotflag = 1 then the ACF is plotted vs lag (sec or m) (default)
%	         2 then the ACF is plotted vs lag
%	         3 then the ACF is plotted vs lag and vs lag (sec or m)
%  for two-dim:
%         L = [Lx Lt] Maximum lag for which the ACF is plotted.
%         (Not implemented!)
%     plotflag = 1 then contour plot (default)
%                2 then mesh plot
%                3 one plot with R(x,0) and one with R(0,t)
%
%   lintype : specify color and lintype, see PLOT for possibilities.
%
% Example:
%   S1 = demospec; S2 = demospec('dir');
%   R1 = spec2cov(S1); R2 = spec2cov(S2,0,30,2,100,[],2,[]);
%   subplot(211), L1 = 17; covplot(R1,L1,1,'.-')
%   subplot(212), covplot(R2)
% 
% See also  dat2cov, dat2cor

% Tested on: matlab 6.0, 5.3, 5.2, 5.1
% History:
% revised by jr 02.04.2000 
% revised by es 24.01.2000 various small things for two-dim  
% revised by es 30.09.1999 (two-dim), 13.10.1999 (ishold and linetype)
% by pab 11.08.98

if ndims(rn.R)>2
  error('Can only handle one- and two dimensional ACF')
end
if (nargin <3 || isempty(plotflag))
  plotflag=1;
end
if nargin<4||isempty(lintype)
  lintype='b-';
end
if ~isfield(rn,'norm')
  rn.norm=0;
end

n=length(rn.R);

if numel(rn.R)==n, % one-dim ---------------------------------
  ih=ishold;

  if n<2, 
    error('The vector must have more than 2 elements!')
  end
  
  if (nargin <2 || isempty(L)),
    L=n-1;
  else
    L=min(L,n-1);
  end
  names=fieldnames(rn);
  ind=find(strcmp(names,'x')+strcmp(names,'t')); %options are 't' and 'x'
  vari=lower(names{ind});
  t=rn.(vari);  %t is either R.t or R.x
  dT=t(2)-t(1);
  if strcmpi(vari,'t')
    unit=' (sec)';
  else
    unit=' (m)';
  end
  
  if isfield(rn,'stdev') && ~isempty(rn.stdev),
    Stdev=rn.stdev(1:(L+1));
    %size(Stdev)
  else
    Stdev=[];
  end
  
  
  if (dT==1 ), % sampling frequency 1 Hz
    plotflag=1; % sufficient to make one plot
  end  
  
  if (plotflag == 3), subplot(211),  end
  
  r = rn.R(1:(L+1));
  tau =  (0:L) ; % # of lags
  mx=max(abs(r));
  if (plotflag == 2) ||(plotflag == 3)
    
    plot(tau,r,lintype) , hold on
    plot([tau(1) tau(L+1)],[0 0],'k-')
    if  ~isempty(Stdev),
      plot(tau,2*Stdev,'r--')
      plot(tau,-2*Stdev,'r--')
    end ,if ~ih, hold off, end
    if dT==1,
      xlabel(['Lag Lag',unit])
    else
      xlabel('Lag')
    end
    ylabel('ACF')
    grid on
    if ih, a=axis; else a=zeros(1,4); end
    axis([0 max(tau(L+1),a(2)) min(-1.01*mx,a(3)) max(1.01*mx,a(4)) ])
    if rn.norm,
      title('Auto Correlation Function (ACF)' )
    else
      title('Auto Covariance Function (ACF)' )
    end
  end
  
  if (plotflag == 3) 
    subplot(212)
  end

  if (plotflag == 1) ||(plotflag == 3)
    tau = tau*dT ;
    plot(tau,r,lintype),hold on
    plot([tau(1) tau(L+1)],[0 0],'k-')
    if  ~isempty(Stdev),	
      plot(tau,2*Stdev,'r--')
      plot(tau,-2*Stdev,'r--')
    end, if ~ih,hold off,end
    xlabel(['Lag',unit])
    ylabel('ACF')
    grid on
    if ih, a=axis; else a=zeros(1,4); end
    axis([0 max(tau(L+1),a(2)) min(-1.01*mx,a(3)) max(1.01*mx,a(4)) ])
    if rn.norm,
      title('Auto Correlation Function (ACF)' )
    else
      title('Auto Covariance Function (ACF)' )
    end
  end
  %subplot(111)
  if ih, hold on, end
else % two-dim ----------------------------------------------------

  [nx,nt]=size(rn.R);
  if (nargin <2 || isempty(L)),
    L=[nx-1, nt-1];
  elseif length(L)==1
    L=[min(L,nx-1) min(L,nt-1)];
  else
    L=[min(L(1),nx-1) min(L(2),nt-1)];
  end
  
  names=fieldnames(rn);
  ind=find(strcmp(names,'x')+strcmp(names,'t')+strcmpi(names,'y'));
  %options are 'x' and 't' and 'y'
  vari1=names{ind(1)}; % always 'x' or 'y' 
  vari2=names{ind(2)}; % always 'y' or 't'
  if strcmpi(vari2,'t')
    unit=' [sec]';
  else
    unit=' [m]';
  end
  x=eval(['rn.',vari1]);
  t=eval(['rn.',vari2]);
  if strcmp(vari2,'y')
    rn.R=rn.R'; % otherwise sizes does not match
  end
  if plotflag == 1
    if strcmp(rn.type,'polar')
      if rn.y(1)>=0 % exploit symmetry
        c=contour([-rn.y(end:-1:2);rn.y],rn.x,[rn.R(:,end:-1:2) rn.R]);
      else
        c=contours(rn.y,rn.x,rn.R);
      end
      limit = size(c,2);
      ix = 1;
      while(ix < limit)
        z_level(ix) = c(1,ix);
        npoints = c(2,ix);
        nexti = ix+npoints+1;
        c(:,ix)=NaN;
        ix = nexti;
      end
      polar(c(1,:),c(2,:),'b');
    else     
      contour(x,t,rn.R',lintype)
      xlabel([vari1 ' [m]'])
      ylabel([vari2 unit])
    end
    if rn.norm,
      title('Auto Correlation Function (ACF)' )
    else
      title('Auto Covariance Function (ACF)' )
    end
  elseif plotflag==2
    surf(x,t,rn.R')
    shading interp
    if rn.norm,
      title('Auto Correlation Function (ACF)' )
    else
      title('Auto Covariance Function (ACF)' )
    end
    xlabel([vari1 ' [m]'])
    ylabel([vari2 unit])
  else
    subplot(211)
    plot(x,rn.R(:,t==0),lintype)
    ih=ishold;
    hold on
    plot([x(1) x(end)],[0 0],':')
    if ~ih, hold off, end
    if rn.norm,
      title('Auto Correlation Function (ACF)' )
    else
      title('Auto Covariance Function (ACF)' )
    end
    xlabel([vari1 ' (m)'])
    subplot(212)
    plot(t,rn.R(x==0,:),lintype)
    ih=ishold;
    hold on
    plot([t(1) t(end)],[0 0],':')
    if ~ih, hold off, end
    if rn.norm,
      title('Auto Correlation Function (ACF)' )
    else
      title('Auto Covariance Function (ACF)' )
    end
    xlabel([vari2 unit])
  end
end
wafostamp()
