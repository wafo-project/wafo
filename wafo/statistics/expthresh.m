function [u,scale,pvalue] = expthresh(data,varargin)
%EXPTHRESH Automatic Threshold selection for exponential data
%
% CALL  [u,shape,scale,pvalue] = expthresh(data,options)
%
% u       = selected threshold
%   scale = WDATA objects with estimated GPD parameters as function of
%           threshold.
% data    = vector of data.
% options = options structure defining the range of thresholds.
%        .method : Method used in the fitting
%        .u    : (default linspace(umin,umax,Nu)
%        .Nmin : Minimum number of extremes to include. (Default Nmin = 10).
%        .umin : Minimum threshold (default min(data))
%        .umax : Maximum threshold (default max(data(end-Nmin)))
%        .Nu   : number of threshold values (default min(length(data),20))
%        .alpha: Confidence coefficient (default 0.05)
%
% EXPTHRESH estimate exponential model parameters over a range of thresholds and
% automatically selects the most appropriate one.
% The purpose is to determine the threshold where the upper tail of the data 
% can be approximated with the exponential distribution. The 
% exponential is appropriate for the tail, if the scale- parameter is constant.
%
% Example
%  opt = expthresh('defaults');
%  opt.Nu = 20;
%  R = rndexp(1,100,1);
%  [u,scale,pvalue] = expthresh(R,opt); 
%  plot(scale); figure(gcf+1);
%  plot(pvalue);
%
%  close all;
%
% See also fitexp, reslife, disprsnidx

%
%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU Lesser General Public License as published by
%     the Free Software Foundation; either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU Lesser General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.


%disp('genparthresh')

  options = struct('method','ML','u',[],'umin',[],'umax',[],'Nu',10,'Nmin',10,'alpha',0.05,'plotflag',0);

  if ischar(data) && strcmpi(data,'defaults')
    u = options;
    return
  end
  options = parseoptions(options,varargin{:});

  sd = sort(data);
  n = length(data);
  if isempty(options.u)
  if isempty(options.Nmin)
     options.Nmin=10;
  else
    options.Nmin = max(options.Nmin,0);
  end
  if options.Nmin>n/2
    warning('WAFO:EXPTHRESH','Nmin too large!')
    options.Nmin = min(options.Nmin,ceil(n/2));
  end
  try
    sdmax = sd(max(n-options.Nmin,5));
  catch
    sdmax = sd(max(n-1,1));
  end
  if isempty(options.umax)
    options.umax = sdmax;
  else
    options.umax = min(options.umax,sdmax);
  end
  xMedian = sd(max(floor(n/2),1));
  if isempty(options.umin)
    options.umin = xMedian;
  else
    options.umin = max(options.umin,sd(1));
  end


  if isempty(options.Nu)
    options.Nu = min(n-options.Nmin,20);
  end

    u1 = linspace(options.umin,options.umax,options.Nu).';
  else
    u1 = options.u(:);
  end
  Nu = length(u1);

  nan1 = nan;
  pvalue1 = nan1(ones(Nu,1));
  phat1 = nan1(ones(Nu,1));
  lo = phat1;
  up = phat1;
  modPvalue = pvalue1;
  okPhat = phat1;

  num = ones(Nu,1);
  method = options.method;

  %xmax = sd(n);
  %returnPrbs = fliplr([1/(n*200) logspace(-log10(n*20),log10(.5))]);
  %RL = invgenpar(returnPrbs/n,0.1,2,2,'lowertail',false);
  %semilogx(1./returnPrbs,RL,'r.'), hold on

  timeSpan = n;
  Tr = 20*timeSpan;
  %
  useModPvalue = 0;

  tmp = zeros(Nu,length(Tr));
  tmplo = tmp;
  tmpup = tmp;




  for ix=Nu:-1:1
    
    [phat] = fitexp(sd(sd>u1(ix))-u1(ix),'method',method,'alpha',options.alpha);
     phat1(ix) = phat.params(1);
     lo(ix) = phat.lowerbound;
     up(ix) = phat.upperbound;   
     pvalue1(ix) = phat.pvalue;
     okPhat(ix) = max(lo(ix:end))<phat1(ix) & phat1(ix)<min(up(ix:end));
     
     
     num(ix) = numel(phat.data);%sum(data>u1(ix));
    if useModPvalue
        modPvalue(ix) = pvalue1(ix);
        returnPrbs = (1:num(ix)).'/(num(ix)+1);
        [RL,RLlo,RLup] = invexp(returnPrbs,phat);
        %excess = max(max(phat.data-RLlo,0),max(RLup-phat.data,0));
        %if any(excess>0)
        excess = phat.data<RLlo|RLup<phat.data;
        if any(excess)
            modPvalue(ix)=exp(log(pvalue1(ix))-sum(excess.*returnPrbs>0));
        end
    end
     
     lambda = num(ix)/timeSpan; % # Y>threshold per timeunit
     returnPrbs = 1./(lambda*Tr);
     [RL,RLlo,RLup] = invexp(returnPrbs,phat,'lowertail',false);
     RL = RL+u1(ix);
     RLlo = RLlo + u1(ix);
     RLup = RLup + u1(ix);
     
    
     tmp(ix,:) = RL;
     tmplo(ix,:) = RLlo;
     tmpup(ix,:) = RLup;
    
    
    
  %    if okPhat(ix) & false
  %     F = edf(phat.data+u1(ix),'wdata',true);
  %     %F.data = (F.data*Nie +0.5)/(Nie+1);
  %     semilogx(Tr,RL,'g-.',Tr,RLlo,'g-.',Tr,RLup,'r-.',1./((1-F.data)*lambda),F.args,'r.');%,...
  %     %semilogx(1./returnPrbs,returnLevels,'r');%,1./returnPrbs,qlo,'b',1./returnPrbs,qup,'b'), hold on
  %     hold('on')
  %     
  %      
  %       shg
  %      %u(ix)
  %      pause(1)
  %      semilogx(1./((1-F.data)*lambda),F.args,'b.');%
  %    end
  end

  if useModPvalue
      pvalue1 = modPvalue;
  end


  dim = 1;
  ixOK = find(okPhat & pvalue1>options.alpha);
  if ~isempty(ixOK)
      if numel(ixOK)>2
        tmp1 = tmp(ixOK,:);
        mRL = mean(tmp1,dim);
        ix12 = findcross(tmp1,mRL);
        ix = ix12(1);
      else
        ix = 1;
      end
      ix = ixOK(ix);
  else
      ix0 = find(okPhat,1);
      if isempty(ix0)
          ix0=1;
      end
      ix = ix0-1+find(pvalue1(ix0:end)>options.alpha,1);
      
      if isempty(ix)
          [maxpvalue,ix] = max(pvalue1(ix0:end));
      end
  end
  u = u1(ix);



  if nargout>1
    p = 1-options.alpha;

    titleTxt2 = sprintf('EXP scale parameter with %d%s CI',100*p,'%');
    titleTxt3 = sprintf('P-value of fitted EXP');


    scale = createwdata('data',phat1(:,1),'args',u1,...
        'dataCI',[lo(:),up(:)],'title',titleTxt2,'labels',{'Threshold','Scale'},...
        'workspace',options,'note',titleTxt2);

    pvalue = createwdata('data',pvalue1,'args',u1,...
        'title',titleTxt3,'labels',{'Threshold','P-value'},...
        'workspace',options,'note',titleTxt3);

    
    scale = wdata(scale);
    pvalue = wdata(pvalue);
    if options.plotflag>0
        subplot(2,1,1)
        plot(pvalue,'.')
        subplot(2,1,2)
        plot(scale,options.plotflag,'.')
    end

  end
end  % main


