function [ind, indg]=findoutliers(xx,zcrit,dcrit,ddcrit,plotflag) 
% FINDOUTLIERS Finds the indices to spurious points in a timeseries
%
% CALL:  [inds, indg] = findoutliers(xn,zcrit,dcrit,ddcrit);
%
%     inds = indices to spurious points.
%     indg = indices to the rest of the points.
% 
%      xn  = two column data matrix with sampled times and values.
%
%    zcrit = critical distance between consecutive points.  
%                 (Default=0)
%    dcrit = critical distance of Dx used for determination of spurious
%            points.  (Default=1.5 standard deviation of xn)            
%   ddcrit = critical distance of DDx used for determination of spurious
%            points.  (Default=1.5 standard deviation of xn)
% plotflag = 0 no plotting (default)
%            1 plot the result
%
%  Consecutive points less than zcrit apart  are considered as spurious.
%  The point immediately after and before are also removed. Jumps greater than
%  dcrit in Dxn and greater than ddcrit in D^2xn are also considered as spurious.
%  (All distances to be interpreted in the vertical direction.)  
%  Another good choice for dcrit and ddcrit are:
%  
%        dcrit = 5*dT  and ddcrit = 9.81/2*dT^2
%
% where dT is the timestep between points.
%
% Example
% xx = load('sea.dat');
% dt = diff(xx(1:2,1));
% dcrit = 5*dt;
% ddcrit = 9.81/2*dt*dt;
% zcrit = 0;
% [inds, indg] = findoutliers(xx,zcrit,dcrit,ddcrit);
% waveplot(xx,'-',xx(inds,:),1,1,1);
%
% assert(length(inds), 1152);
% close all;
%
% See also  waveplot, reconstruct

% Tested on: Matlab 5.3, 5.2 , 5.1
% History:           
% last modified by Per A. Brodtkorb 
% 29.03.99, new input arguments
% 01.10.98 checks  input and accepts missing values NaN's
% 25.09.98, 13.08-98


if nargin <5||isempty(plotflag)
  plotflag=0;
end
% finding outliers
findjumpsDx=1; % find jumps in Dx
%      two point spikes and Spikes dcrit above/under the
%       previous and the following point are spurios.
findSpikes=0; % find spikes
findDspikes=0; % find double (two point) spikes
findjumpsD2x=1; % find jumps in D^2x
findNaN=1; % find missing values

xn=xx;
[n m]= size(xn);

if n<m
 b=m;m=n;n=b; 
 xn=xn';
end

if n<2, 
  error('The vector must have more than 2 elements!');
end

switch m
 case 1, % OK dimension
 case 2, xn=xn(:,2);
 otherwise, 
   error('Wrong dimension of input! dim must be 2xN, 1xN, Nx2 or Nx1 ');
end

ind=[];indg=[];
indmiss=isnan(xn);


if nargin<2||isempty(zcrit),
  zcrit=0;
end
if nargin<3||isempty(dcrit),
  dcrit=1.5*std(xn(~indmiss));
  disp(['dcrit is set to ' num2str(dcrit)]);
end
if nargin<4||isempty(ddcrit),
  ddcrit=1.5*std(xn(~indmiss));
  disp(sprintf('ddcrit is set to %s', num2str(ddcrit)));
end
if findNaN, 
  ind=find(indmiss);
  disp(sprintf('Found %d missing points', length(ind)));
  xn(indmiss)=0;%set NaN's to zero 
end

dxn=diff(xn);
ddxn=diff(dxn);

if  findSpikes , % finding spurious spikes
  tmp=find(((dxn(1:(end-1))>dcrit).*(dxn(2:end))<-dcrit) |...
      ((dxn(1:(end-1))<-dcrit).*(dxn(2:end))>dcrit) )+1;
  %if ~isempty(tmp),
    disp(sprintf('Found %d spurious spikes', length(tmp)));
  %end
  ind=[ind;tmp];
end

if findDspikes ,% finding spurious double (two point) spikes  
  tmp= find(((dxn(1:(end-2))>dcrit).*(dxn(3:end))<-dcrit) |...
      ((dxn(1:(end-2))<-dcrit).*(dxn(3:end))>dcrit) )+1;
  %if ~isempty(tmp),
    disp(sprintf('Found %d spurious two point (double) spikes', length(tmp)));
  %end
  ind=[ind;tmp;tmp+1];%removing both points
end

if findjumpsDx ,% finding spurious jumps  in Dx
  % finding spurious positive jumps  
  tmp= find(dxn>dcrit);
  %if ~isempty(tmp),
    disp(sprintf('Found %d spurious positive jumps of Dx', length(tmp)));
  %end
  ind=[ind;tmp+1]; %removing the point after the jump 

  % finding spurious negative jumps  
  tmp= find(dxn<-dcrit);
  %if ~isempty(tmp),
    disp(sprintf('Found %d spurious negative jumps of Dx', length(tmp)));
  %end
  ind=[ind;tmp];% tmp+1]; %removing the point before the jump
end

if findjumpsD2x ,% finding spurious jumps in D^2x  
  % finding spurious positive jumps  
  tmp= find(ddxn>ddcrit)+1;
  %if ~isempty(tmp),
    disp(sprintf('Found %d spurious positive jumps of D^2x',length(tmp)));
  %end
  ind=[ind;tmp];%tmp+1];%tmp-2]; removing the jump

  % finding spurious negative jumps  
  tmp= find(ddxn<-ddcrit)+1;
  %if ~isempty(tmp),
    disp(sprintf('Found %d spurious negative jumps of D^2x',length(tmp)));
  %end
  ind=[ind;tmp];%tmp+1];% tmp-2];removing the jump
end

if zcrit>=0
  % finding consecutive values less than zcrit apart.
  indzeros=(abs(dxn)<=zcrit);
  indz=find(indzeros)+1;
  %if ~isempty(indz),
    if zcrit==0,
      disp(sprintf('Found %d consecutive equal values',length(indz)));
    else
      disp(sprintf('Found %d consecutive values',length(indz)));
      disp(sprintf('less than %s apart', num2str(zcrit)));
    end
  %end
  
  %finding the beginning and end of consecutive equal values
  indtr=find((diff(indzeros)))+1;

  %indices to consecutive equal points
  if 1, % removing the point before + all equal points + the point after  
    ind=[ind;(indtr(:)-1);indz;indtr(:);(indtr(:)+1);];
  else % removing all points + the point after 
  %  ind=[ind;indz;indtr(:);(indtr(:)+1)];
  end  
end

if length(ind)>1,
  ind=sort(ind);
  ind(diff(ind)==0)=[];%removing indices to points identified several times
end

%if ~isempty(indg),
  % for ix=1:length(indg)
  % ind(indg(ix)==ind)=[];%removing indices to points near NaN's
  % end
%end
if ~isempty(ind),
  disp(sprintf('Found the total of %d spurious points', length(ind)));
end

if (nargout==2)||plotflag
  indg=(1:n)';
  indg(ind)=[];
end


 if plotflag,
   %xn2=[xn(:,1)   spline(xn(indg,1),xn(indg,2),xn(:,1)) ]; 
   waveplot(xx(indg,:),xx(ind,:),20)
 end
