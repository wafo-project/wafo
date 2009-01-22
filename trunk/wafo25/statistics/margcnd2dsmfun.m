function [Av , Bv , Cv ]=margcnd2dsmfun(phat,h2,CSMA,lin_extrap)
%MARGCND2DSMFUN Smooths the MARGCND2D distribution parameters. 
%
% CALL:  [A, B, C] = margcnd2dsmfun(phat,x2,csm,lin)
%
%   A, B, C = scale shape and location parameter, respectively
%             evaluated at x2
%   phat    = parameter structure array (see fitmargcnd2d)
%   x2      = evaluation points (default phat.x{1}(:,1)) 
%   csm     = [csma csmb csmc], smoothing parameter vector which defines
%             the smoothing of parameter A,B and C, respectively 
%             0 -> LS straight line
%             1 -> cubic spline interpolant (default [1 1 1])
%   lin     = [lina linb linc], vector defining the extrapolation of
%             the parameters A,B and C,  respectively (default [1 1 1])
%             0 No linear extrapolation outside the range of data
%             1 Linear extrapolation outside the range of data 
%
%   The size of A B and C is the size of the argument  X2.
%   
% See also   fitmargcnd2d, margcnd2dsmfun2, smooth


% Tested on matlab 5.x
% history
% revised pab July2004
%  -made useSig work correctl
% revised pab 03.12.2000
% - added truncated weibull and truncated rayleigh
% revised pab 12.11.2000
%  - added pdfgengam option
%  Per A. Brodtkorb 28.10.98

PVH   = phat.x{1};
ind   = find(~isnan(sum(PVH(:,:),2)));% & PVH(:,1)<2);
PVH   = PVH(ind,:);
CDIST = phat.dist{1};  
if nargin<2 ||isempty(h2),
  h2=PVH(:,1);
end

CSMB=1;
CSMC=1;
if nargin<3||isempty(CSMA)
  CSMA = 1;
elseif length(CSMA)==2
  CSMB=CSMA(2);
  CSMA=CSMA(1);
elseif length(CSMA)>=2
  CSMB=CSMA(2);
  CSMC=CSMA(3);
  CSMA=CSMA(1);
end
linB=1;
linC=1;
if (nargin<4) || isempty(lin_extrap),
    linA=1;
elseif length(lin_extrap)==2,
  linA=lin_extrap(1);
  linB=lin_extrap(2);
elseif length(lin_extrap)>=2,
  linA=lin_extrap(1);
  linB=lin_extrap(2);
  linC=lin_extrap(3);
end

%smooth extrapolates linearly outside the range of PVH
useSIG = (1 && isfield(phat,'CI')) ; %
useSTATS = (useSIG && isfield(phat,'stats1') && ~isempty(phat.stats1));
if useSTATS
  Npts =  phat.stats1{end}(ind);
end
if useSIG
  useSIG = any(~isnan(phat.CI{1}(:)));
end
K5=1;
if all(PVH(:,1)>=0)
  da = (PVH(:,1)).^K5.*ones(size(PVH(:,1)));
else
  da = ones(size(PVH(:,1)));
end
if useSIG, 
  CIa = phat.CI{1}(ind,1:2);
  if useSTATS
    NptsA = Npts;
  end
  k1 =  find(PVH(:,2)<CIa(:,1)| CIa(:,2)< PVH(:,2));
  if any(k1)
    CIa(k1,:) = [];
    PVH(k1,:)=[];
    ind(k1) = [];
    da(k1)  = [];
    if useSTATS
      NptsA(k1) = [];
    end
  end
  tmp = (diff(CIa,1,2)/4).^2; % approximate 2 standard
                                           % deviations
   if ~isempty(tmp)
      k0 = find(~isnan(tmp));
      maxStd = 10;
      if any(k0)
	da(k0) = tmp(k0);
	maxStd = max(maxStd,10*max(tmp(k0)));
      end
      k1 = find(isnan(tmp));
      if any(k1)
	da(k1) = maxStd;
      end
      %size(PVH),size(da)
      
       da = da+[0; da(1:end-1)]+...
	    [da(2:end);da(end)]+...
       PVH(:,1).*[0; abs(diff(PVH(:,2)))];
       if useSTATS
	 da = da./sqrt(NptsA);
       end
      da = PVH(:,1).*da./(PVH(:,2)+eps);
   end				     
end

Av=smooth(PVH(:,1),PVH(:,2),CSMA,h2,linA,da);
switch lower(CDIST(1:2)),
  case {'ra','tr','tg','gu','ga','we','tw','gg'}, 
    if (any(Av<0)),
      da2 = abs(log(PVH(:,2)) - log(PVH(:,2)+ da));
      Av=exp(smooth(PVH(:,1),log(PVH(:,2)),CSMA,h2,linA,da2));
    end
end
Bv=[];	  
if ~strcmpi(CDIST(1:2),'ra')
  if all(PVH(:,1)>=0)
    db = PVH(:,1).^K5.*ones(size(PVH(:,1)));
  else
    db = ones(size(PVH(:,1)));
  end
  if useSTATS
    NptsB = Npts;
  end
  if useSIG 
    CIb = phat.CI{1}(ind,3:4);
    k1 =  find(PVH(:,3)<CIb(:,1)| CIb(:,2)< PVH(:,3));
    if any(k1)
      CIb(k1,:) = [];
      PVH(k1,:)=[];
      ind(k1) = [];
      db(k1) = [];
      if useSTATS
	NptsB(k1) = [];
      end
    end
    
    % approximate 2 standard deviations
    tmp=(diff(CIb,1,2)/4).^2; 
    if ~isempty(tmp)
      k0 = find(~isnan(tmp));
      maxStd = 10;
      if any(k0)
        db(k0) = tmp(k0);
        maxStd = max(maxStd,10*max(tmp(k0)));
      end
      k1 = find(isnan(tmp));
      if any(k1)
        db(k1) = maxStd;
      end
      db = (db+[0; db(1:end-1)]+...
        [db(2:end);db(end)])/3+...
        PVH(:,1).*[0; abs(diff(PVH(:,3)))];
      if useSTATS
        db = db./sqrt(NptsB);
      end
      db = PVH(:,1).*db./(PVH(:,3)+eps);
      %db = db./(PVH(:,3)+eps);
    end
  end
  
  Bv=smooth(PVH(:,1),PVH(:,3),CSMB,h2,linB,db);
  switch CDIST(1:2) ,
    case {'lo','ga','we','tw','gg'},  
      if  (any(Bv<0)),
	db2 = abs(log(PVH(:,3)) - log(PVH(:,3)+ db));
	Bv=exp(smooth(PVH(:,1),log(PVH(:,3)),CSMB,h2,linB,db2));
      end
  end 
end

Cv=0;
if (((size(PVH,2)>=4) && (~strcmpi(CDIST(1:2),'ra'))) || ...
    ((size(PVH,2)>=3) && (strcmpi(CDIST(1:2),'ra')))),
  Cv=smooth(PVH(:,1),PVH(:,4),CSMC,h2,linC);
  
  if all(PVH(:,1)>=0)
    dc = PVH(:,1).^K5.*ones(size(PVH(:,1)));
  else
    dc = ones(size(PVH(:,1)));
  end

  if useSIG
    CIc = phat.CI{1}(ind,5:6);
    if useSTATS
      NptsC = Npts;
    end
    k1 =  find(PVH(:,4)<CIc(:,1)| CIc(:,2)< PVH(:,4));
    if any(k1)
      CIc(k1,:) = [];
      PVH(k1,:)=[];
      %ind(k1) = [];
      dc(k1) = [];
      NptsC(k1) = [];
    end
    tmp = (diff(CIc,1,2)/4).^2;
    % approximate 2 standard deviations
    if ~isempty(tmp)
      k0 = find(~isnan(tmp));
      maxStd = 10;
      if any(k0)
        dc(k0) = tmp(k0);
        maxStd = max(maxStd,10*max(tmp(k0)));
      end
      k1 = find(isnan(tmp));
      if any(k1)
        dc(k1) = maxStd;
      end
      dc = dc+[0; dc(1:end-1)]+...
        [dc(2:end);dc(end)];%+...
      %PVH(:,1).*[0; abs(diff(PVH(:,4)))];;
       if useSTATS
         dc = dc./sqrt(NptsC);
       end
       dc = PVH(:,1).*dc./(PVH(:,4)+eps);
    end
  end
  switch lower(CDIST(1:2)) ,
    case {'ra','tg','lo','ga','we'},  
      if  (any(Cv<0)),
        Cv(Cv<0)=0;
      end
    case 'gg',
       if  (any(Cv<0)),
         Cv=exp(smooth(PVH(:,1),log(PVH(:,4)),CSMC,h2,linC,dc));
      end 
  end 
end
return


