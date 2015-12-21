function Nf1 = waveplot(x,varargin)
%WAVEPLOT Plots the surface elevation of timeseries.
%
% CALL:  waveplot(x1,x2,Nsub,Nf,Sm0,v_fact,sym1,sym2)
%
%     x1, x2 = two-column timeseries 
%              first column sampling times [sec]
%              second column surface elevation [m]
%              x2 is by default zero-separated troughs and crests.
%       Nsub = Number of subplots in each figure. By default 
%              Nsub is such that there are about 20 mean down 
%              crossing waves in each subplot. If Nf is not 
%              given and Nsub is larger than 6 then Nsub is	       
%              changed to Nsub=min(6,ceil(Nsub/Nf))
%        Nf  = Number of figures. By default Nf=ceil(Nsub/6). 
%       Sm0  = standard deviation of x1. 
%     v_fact = how large in stdev the vertical scale should be (default 3)
% sym1, sym2 = plot symbol and color for x1 and x2, respectively 
%              (see PLOT)  (default 'k.' and 'k+') 
%
%  Note: - sym1 and sym2 can be given anywhere after x1.
%          if omitted default values are used.
%        - x2 can  be omitted, but if given it must appear 
%          before the scalars Nsub,Nf,Sm0 and v_fact. 
%        - if [] is given for any of the scalars default values are used.
%
% Example: 
% % Plot x1 with red lines and mark troughs and crests with 
% % blue circles.
%   x = load('sea.dat');
%   x1 = x(1:150,:);
%   waveplot(x1,'r-','bo')
%
% See also  dat2tc, plot

%Tested on: Matlab 6.0, 5.3, 5.2, 5.1
% History:
% revised jr 02.04.2001
%  - added example, updated info.
% revised pab 11.10.2000
%  - bug fix: disabled subplot when Nsub=1 => may use waveplot to
%             overplot other figures 
% revised pab 01.02.2000
%  - fixed a bug in hour scale
%  - temporarily disabled this function
% revised pab 03.12.1999
%  - changed input to varargin, removed Nw
%  - added hour and minute scale for horizontal axis i.e. dT
% last modified Per A. Brodtkorb 01.10.98 
%    to accept missing values NaN's
% revised pab 15.08.98
% 
xn=x;

[n m]= size(xn);
if n<m
 b=m;m=n;n=b; 
 xn=xn';
end

if n<2, 
  error('The vector must have more than 2 elements!')
end

istime=1;

switch m
 case 1, xn=[ (1:n)' xn(:)];istime=0;
 case 2, % dimension OK!
 otherwise, error('Wrong dimension of input! dim must be 2xN, 1xN, Nx2 or Nx1 ')          
end

[TP1,Nsub,Nf,Sm0,sym1,sym2, v_fact] = wavechk(varargin,xn);

Ns = floor(n/(Nf*Nsub));
ind = 1:Ns;
if all(xn(:,2)>=0)
  v_scale = [0 2*Sm0]*v_fact;
else
   v_scale = [-Sm0 Sm0]*v_fact;
end
if istime
  XlblTxt='Time (sec) ';dT=1;
  if 1, % disable other scalings pab 01.02.2000
    if abs(xn(ind(1),1)-xn(ind(Ns),1))>18000, % more than 5 hours
      dT=1/(60*60);
      XlblTxt='Time (hours) ';
    elseif abs(xn(ind(1),1)-xn(ind(Ns),1))>300,% more than 5 minutes
      dT=1/60;
      XlblTxt='Time (minutes) ';
    end
  end   
else
  dT=1;
  XlblTxt='Step ';
end 

indmiss=isnan(xn(:,2)); % indices to missing points
if max(abs(xn(~indmiss,2)))>5*Sm0,
    XlblTxt=[ XlblTxt '(Spurious data since max > 5 std.)'];
end
start=gcf-1; %  start at current figure
hstate = ishold;

for iz=1:Nf
  if Nf>1,figure(start+iz);end
  for ix=1:Nsub,
    if Nsub>1,subplot(Nsub,1,ix),end
    if hstate, hold on; end
    h_scale=[xn(ind(1),1) xn(ind(Ns),1)]*dT;
    plot(xn(ind,1)*dT,xn(ind,2:end),sym1) , hold on
    plot(TP1(:,1)*dT,TP1(:,2:end),sym2) 
    plot(h_scale,[0 0],'k')
    axis([h_scale v_scale])
  
    for iy=[-2 2],
      plot(h_scale,iy*Sm0*[1 1],'--')
      switch sign(iy)
        case -1, sign1='-';
        case  0, sign1=' ';
        case  1, sign1='+';
      end
      figtext(0.99,iy*Sm0,['- ' sign1 num2str(abs(iy)) ],'norm','data')
    
    end, hold off
  
    ind = ind + Ns; 
  end

  xlabel(XlblTxt)

  if Nsub>1,subplot(Nsub,1,1),end

  title('Surface elevation from mean water level (MWL).')

  if Nsub>1,subplot(Nsub,1,floor((Nsub+1)/2)),end
  if (Sm0 >1.1)||(Sm0<0.9), % surface elevation is probably not standardized 
    ylabel('Distance from MWL.(m)')
  end
  set(gca,'DefaultTextRotation',90)
  figtext(1.05,0,'Standard Deviation','norm','data','center')
  set(gca,'DefaultTextRotation',0)
 wafostamp 
end

if nargout>0,
  Nf1=Nf;
end

%subplot(111)
return


function [x2,Nsub,Nf,Sm0,sym1,sym2, vfact] = wavechk(P,x1)
%WAVECHK Helper function for waveplot.
%
% CALL  [x2, Nsub,Nf,Sm0,sym1,sym2, vfact]=wavechk(P,x1) 
%
%   P = the cell array P of input arguments (between 0 and 7 elements)
%  x1 = must be a two column vector.


Nw=20;
Np=length(P);
strix=zeros(1,Np);
for ix=1:Np, % finding symbol strings 
 strix(ix)=ischar(P{ix});
end
sym1='k.'; % Black dots is default
sym2='k+'; % Black plus is default
k=find(strix);
if any(k)
  sym1=P{k(1)};
  Np=Np-1;
  if length(k)>1
    sym2=P{k(2)};
    Np=Np-1;
  end
end

P={P{find(~strix)}};

indmiss=isnan(x1(:,2)); % indices to missing points
if (Np>0) && (numel(P{1})>1)
  x2=P{1};
  P={P{2:Np}};
  Np=Np-1;
else
  x2=dat2tc(x1(~indmiss,:),0,'tw'); % Finding mean separated turning points 
  %TP1=data2tp(xn,0.2,'mw'); % Finding rf waves
  %TP1=dat2tp(xn,0.2);
end

if (Np >= 1) && ~isempty(P{1})
% waveplot(x1,x2,Nsub)
  Nsub=P{1};
else
  Nsub=floor(length(x2)/(2*Nw))+1; % about Nw mdc waves in each plot
end

if (Np >= 2) && ~isempty(P{2})
  % waveplot(x1,x2,Nsub,Nf)
  Nf=P{2};
else
  Nf=ceil(Nsub/6); 
  Nsub=min(6,ceil(Nsub/Nf));
end
if (Np >= 3) && ~isempty(P{3})
  % waveplot(x1,x2,Nsub,Nf,Sm0)
  Sm0=P{3};
else
  Sm0=std(x1(~indmiss,2));
end
if (Np >= 4) && ~isempty(P{4})
   % waveplot(x1,x2,Nsub,Nf,Sm0,v_fact)
 vfact=P{4};
else
  vfact=3; % 3 std default
end
