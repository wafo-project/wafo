function spwaveplot(xx,varargin)
%SPWAVEPLOT Plots specified waves from a timeseries
%
% CALL:  spwaveplot(x1,x2,w_ind,tz_ind,sym1,sym2)
%
%      x1,x2 =  two-column timeseries 
%               first column sampling times [sec]
%               second column surface elevation [m]
%
%      w_ind = vector of indices to waves we want to plot, i.e., wave numbers. 
%     tz_ind = vector of indices to the beginning, middle and end of 
%              defining wave, i.e. for zero-downcrossing waves, indices to 
%              zerocrossings (default trough2trough wave)
% sym1, sym2 = plot symbol and color for x1 and x2, respectively 
%               (see PLOT)  (default 'k-' and 'k+') 
%
%  Note: - sym1 and sym2 can be given anywhere after x1.
%           if omitted default values are used.
%         - x2 can  be omitted, but if given it must appear 
%           before the vectors w_ind, tz_ind
%
% Example: Plot waves nr. 6,7,8 and waves nr. 12,13,...,17
%  x = load('sea.dat'); x1 = x(1:500,:);
%  spwaveplot(x1,[6:8 12:17])
%
% See also  waveplot, dat2tc

% Tested on: Matlab 6.0, 5.3, 5.2, 5.1
% History:
% revised jr 02.04.2001
%  - updated example, documentation
% revised pab 24.05.2000
%  - improved order of input arguments
%  - improved w_ind and plotting
% revised pab 03.12.1999
%  - updated documentation
% by pab 11.08.99

[xx2,Nsub,Nfig,Nwp,w_ind,tz_ind,sym1,sym2] = wavechk(varargin,xx);


for iy=1:Nfig
  for ix=1:Nsub,
    subplot(Nsub,1,mod(ix-1,Nsub)+1)
    ind= tz_ind(2*w_ind(ix)-1):tz_ind(2*w_ind(ix)+2*Nwp(ix)); % indices to wave
    plot(xx(ind,1),xx(ind,2),sym1);hold on %plot wave
    plot(xx([ind(1) ind(end)],1),[0 0]),  % plot zero level
    
    if ~isempty(xx2)
      ax = axis;
      axis([xx([ind(1) ind(end)],1)' ax(3:4)]) %make sure the axis does not change
      plot(xx2(:,1),xx2(:,2),sym2) % before we plot xx2 on top
    end
    hold off
    if Nwp(ix)==1,
      ylabel(['Wave ', int2str(w_ind(ix))])
    else
      ylabel(['Wave ', int2str(w_ind(ix)) '-' int2str(w_ind(ix)+Nwp(ix)-1)])
    end 
    
  end % Nsub
  xlabel('Time (sec)');
  wafostamp
  if (mod(ix-1,Nsub)+1==Nsub)&& iy~=Nfig, figure(gcf+1);  end
end % Nfig


function [x2,Nsub,Nfig,Nwp,w_ind,tz_ind,sym1,sym2] = wavechk(P,x1)
%WAVECHK Helper function for spwaveplot.
%
% CALL  [x2, Nsub,Nf,w_ind,tz_ind,sym1,sym2]=wavechk(P,x1) 
%
%   P = the cell array P of input arguments (between 0 and 7 elements)
%  x1 = must be a two column vector.



Np=length(P);
strix=zeros(1,Np);
for ix=1:Np, % finding symbol strings 
 strix(ix)=ischar(P{ix});
end
sym1='k-'; % Black line is default
sym2='k+'; % Black plus is default
k=find(strix);
if any(k)
  sym1=P{k(1)};
  Np=Np-1;
  if length(k)>1
    sym2=P{k(2)};
    Np=Np-1;
  end
  P={P{find(~strix)}};
end


x2=[];

if (Np>0) && (all(size(P{1})>1))
  x2=P{1};
  P={P{2:Np}};
  Np=Np-1;
end

if (Np >= 1) && ~isempty(P{1})
% waveplot(x1,x2,Nsub)
  w_ind=P{1};
else
  error('w_ind must be given')
end

if (Np >= 2) && ~isempty(P{2})
  tz_ind=P{2};
else
  %indmiss=isnan(x1(:,2)); % indices to missing points
  [tc tc_ind,tz_ind]=dat2tc(x1,0,'tw'); % finding trough to trough waves
end

%Nw = length(w_ind);%# waves
dw = find([0;abs(diff(w_ind(:)))>1]);
Nsub =length(dw)+1;
Nwp=zeros(1,Nsub);
if Nsub>1
  Nwp(Nsub)= w_ind(end)-w_ind(dw(end))+1;
  w_ind(dw(end)+1:end)=[];
  for ix=Nsub-2:-2:1
    Nwp(ix+1) = w_ind(dw(ix+1)-1)-w_ind(dw(ix))+1; % # of waves pr subplot
    w_ind(dw(ix)+1:dw(ix+1)-1)=[];
  end
  Nwp(1)= w_ind(dw(1)-1)-w_ind(1)+1;
  w_ind(2:dw(1)-1)=[];
else
  Nwp(1)= w_ind(end)-w_ind(1)+1;
end

Nsub=min(6,Nsub);
Nfig=ceil(Nsub/6); 
Nsub=min(6,ceil(Nsub/Nfig));

