function [tpe,Pout,I,tpe0] = tpextrapolate(tp,N,Pin,plotflag)

%TPEXTRAPOLATE Extrapolates a sequence of turning points.
%
% CALL: tpe = tpextrapolate(tp,N);
%       [tpe,Pout,I] = tpextrapolate(tp,N,Pin,plotflag)
%
%   tp          = A sequence of turning points.         [n,1] / [n,2]
%   N           = Number of blocks to extrapolate.
%   Pin         = Input parameters. (optional)          [struct array]
%    .method    = Method for extrapolating. 
%                'exp' : Exponential distribution  (default)
%                'gpd' : Generalized Pareto Distribution. 
%                        Note: Use 'gpd' only if 'exp' gives bad fit to data.
%    .LCfrac    = Fraction of level crossings (LC) for automatic choice of levels, u_lev.  
%                 The upper and lower levels for extrapolation are chosen so that they are 
%                 crossed LCmax*LCfrac number of times.  LCmax is the maximum of LC.
%                 Default: Pin.LCfrac=1/sqrt(LCmax)   
%    .u_lev     = Lower and upper levels for extrapolation.  (manual choice of levels)
%                   Pin.u_lev=[u_min u_max]
%    .lim       = Limits of signal,  Pin.lim=[minlim maxlim] .
%                 Values outside the range is set to the limit. 
%                 Default: Pin.lim=[], no limits
%
%   plotflag    = 0: Don't plot diagnostic plots, (default)
%                 1: Plot exponential probability plots,
%                 2: Also plot mean excess plots.
%
%   tpe         = Extrapolated turning points           [N*n,1] / [N*n,2]
%   Pout        = Output parameters. (see Pin)                   [struct array]
%   I           = Indeces to the extrapolated points.
%   
% The highest maxima and the lowest minima of the turning points are extrapolated.  
% The result is an N times as long signal consisting of N blocks.  Each block is 
% generated from tp with randomly simulated maxima above u_max, and minima below u_min.
% The mean of the exceedances above u_max and below u_min are estimated from tp.
%
% Example:
%   x = load('sea.dat');
%   tp = dat2tp(x,0.5);
%   tpe = tpextrapolate(tp,1,[],1);
%   clf, plot(tp(:,1),tp(:,2),'b',tpe(:,1),tpe(:,2),'r')
%   [tpe,Pout,I] = tpextrapolate(tp,1,[],2);
%   clf, plot(tpe(:,1),tpe(:,2),'b',tpe(I.min,1),tpe(I.min,2),'g*',tpe(I.max,1),tpe(I.max,2),'g*')
%
% See also  rfmextrapoalte

% Tested  on Matlab  6.5
%
% History:
% Created by PJ (Pär Johannesson) 16-Apr-2003
% Updated by PJ  11-Jun-2003
% Updated by PJ  24-Jun-2003
% Updated by PJ  05-Sep-2003
%   Added output I, and input Pin.lim.
%   Now also handles zero number of exceedances.
%   Added GPD.
% Updated by PJ  29-Oct-2003
%   If N<1, estimate parameters Pout, but set tpe=[].
% Updated by PJ  11-Mar-2005
%   Added "Conservative" extrapolation of load spectrum.


% Check input arguments
ni = nargin;
no = nargout;
error(nargchk(2,4,ni));

if ni<3, Pin=[]; end
if ni<4, plotflag=[]; end

if isempty(plotflag)
    plotflag=0;
end

%
% Default Parameters
%

Pout.method = 'exp';
Pout.LCfrac = []; %0.05;
Pout.u_lev = [];
Pout.lim = [];

% Copy input parameters Pin to Pout
if ~isempty(Pin)
    Fname = fieldnames(Pin);
    for i = 1:length(Fname)
        Pout = setfield(Pout,Fname{i},getfield(Pin,Fname{i}));
    end
end

method = lower(Pout.method);
if ~(strcmp(method,'exp') | strcmp(method,'gpd'))
    error(['Undefined method: ' Pout.method]), 
end

if isempty(Pout.u_lev)
    lc = tp2lc(tp);  % Level crossings
    
    LCmax=max(lc(:,2));
    if isempty(Pout.LCfrac)
        Pout.LCfrac = 1/sqrt(LCmax);
    end
    
    nLCfrac = LCmax*Pout.LCfrac;
    imax = max(find(lc(:,2)>nLCfrac));
    imin = min(find(lc(:,2)>nLCfrac));
    umax = lc(imax,1);
    umin = lc(imin,1);
    Pout.u_lev = [umin umax];
end

u_min = Pout.u_lev(1);
u_max = Pout.u_lev(2);

if tp(1,end)>tp(2,end)
    StartMax=1; StartMin=2;
else
    StartMax=2; StartMin=1;
end

% Estimate excedances

Imax = 2*(find(tp(StartMax:2:end,end)>u_max)-1)+StartMax;
Imin = 2*(find(tp(StartMin:2:end,end)<u_min)-1)+StartMin;

y_max = tp(Imax,end)-u_max;
y_min = -(tp(Imin,end)-u_min);

a_max = mean(y_max);
a_min = mean(y_min);
if strcmp(method,'gpd')
    gpd_max = fitgenpar(y_max);
    gpd_min = fitgenpar(y_min);
    ci_k_max = [gpd_max.lowerbound(1) gpd_max.upperbound(1)];
    ci_k_min = [gpd_min.lowerbound(1) gpd_min.upperbound(1)]

end

if plotflag>0
    subplot(2,2,1),plot(1:length(y_max),y_max,'.'), 
    title(['Exceedances of maxima above u_{max}=' num2str(u_max)]), ylabel('Exceedances of maxima')
    subplot(2,2,2), if length(y_max)>0, plotexp(y_max), end
    subplot(2,2,3),plot(1:length(y_min),y_min,'.'), 
    title(['Exceedances of minima below u_{min}=' num2str(u_min)]), ylabel('Exceedances of minima')
    subplot(2,2,4), if length(y_min)>0, plotexp(y_min), end
    drawnow
    
    Pout.Zmin = y_min;
    Pout.Zmax = y_max;
end

if N>0
    meth=2;  % Method 2
    conservarive = 1;  % 1=Conservative extrapolation
    
    if meth==2
        [yy_max,IImax] = sort(y_max);
        [yy_min,IImin] = sort(y_min);
    end
    tpe = [];
    for k = 1:N
        
        % Simulate independent Exp or GPD exceedances
        
        if length(y_max)>0
            if strcmp(method,'exp')
                yr_max = rndexp(a_max,length(y_max),1);
            else
                yr_max = rndgenpar(gpd_max,length(y_max),1);
            end
        else
            yr_max=[];
        end
        
        if length(y_min)>0
            if strcmp(method,'exp')
                yr_min = rndexp(a_min,length(y_min),1); 
            else
                yr_min = rndgenpar(gpd_min,length(y_min),1);
            end
        else
            yr_min=[];
        end
        
        if meth ==1
            % Method 1
            % Independent ordering of excedances
            tpr = tp;
            tpr(Imax,end) = u_max + yr_max;
            tpr(Imin,end) = u_min - yr_min;
        else
            
            % Method 2
            % Simulate independent Exp
            % Order the sample according to the order of the measurements.
            
            yr_max = sort(yr_max);
            yr_min = sort(yr_min);
            
            tpr = tp;
            tpr(Imax(IImax),end) = u_max+yr_max;
            tpr(Imin(IImin),end) = u_min-yr_min;
        end
        
        if conservarive
            I=find(tpr(Imin,end)>tp(Imin,end));
            tpr(Imin(I),end) = tp(Imin(I),end);
            I=find(tpr(Imax,end)<tp(Imax,end));
            tpr(Imax(I),end) = tp(Imax(I),end);
        end
        
        tpe = [tpe; tpr];
    end
    
    if no>3
        tpe0=tpe;
    end
    
    % Apply limits
    if ~isempty(Pout.lim)
        minlim=Pout.lim(1);  maxlim=Pout.lim(2);
        if isnumeric(maxlim) % Don't apply if maxlim=Inf
            I = find(tpe(:,end)>maxlim);
            if ~isempty(I)
                tpe(I,end)=maxlim; 
            end
        end
        if isnumeric(minlim) % Don't apply if minlim=-Inf
            I = find(tpe(:,end)<minlim);
            if ~isempty(I)
                tpe(I,end)=minlim;
            end
        end
    end
    
    % Make sure that the output is a sequence of turning points
    tpe = rfcfilter(tpe,0,1);
    %tpe = dat2tp(tpe);  
    
else
    tpe = [];
end

% Output parameters
Pout.a_min = a_min;
Pout.a_max = a_max;
Pout.n_min = length(y_min);
Pout.n_max = length(y_max);

if strcmp(method,'gpd')
    Pout.gpd_min = gpd_min;
    Pout.gpd_max = gpd_max;
    Pout.ci_k_min = ci_k_min;
    Pout.ci_k_max = ci_k_max;
end

if N>0
    % Indeces of extrapolated values
    if no>2
        if tpe(1,end)>tpe(2,end)
            StartMax=1; StartMin=2;
        else
            StartMax=2; StartMin=1;
        end
        
        % Estimate excedances
        I=[];
        I.max = 2*(find(tpe(StartMax:2:end,end)>u_max)-1)+StartMax;
        I.min = 2*(find(tpe(StartMin:2:end,end)<u_min)-1)+StartMin;
    end
end

% Diagnostic plots 2
if plotflag>1
    
    m = mean(tp(:,end));
    
    n=500;
    Umin = linspace(min(tp(:,end)),m,n);
    Umax = linspace(m,max(tp(:,end)),n);
    Amin=zeros(1,n); Amax=zeros(1,n);
    dAmin=zeros(1,n); dAmax=zeros(1,n);
    for k=1:n
        Imax = 2*(find(tp(StartMax:2:end,end)>Umax(k))-1)+StartMax;
        Imin = 2*(find(tp(StartMin:2:end,end)<Umin(k))-1)+StartMin;
        
        y_max = tp(Imax,end)-Umax(k);
        if ~isempty(y_max), Amax(k) = mean(y_max); dAmax(k)=Amax(k)/sqrt(length(y_max)); end
        y_min = -(tp(Imin,end)-Umin(k));
        if ~isempty(y_min), Amin(k) = mean(y_min); dAmin(k)=Amin(k)/sqrt(length(y_min)); end
    end
    
    figure
    subplot(2,1,1)
    plot(Umax,Amax,'r',Umax,Amax-2*dAmax,'b:',Umax,Amax+2*dAmax,'b:'), grid
    title('Maxima - Choose a level where the estimate is stable'), 
    ylabel('Estimated mean exceedance, m'), xlabel('Upper threshold level, u_{max}')
    subplot(2,1,2)
    plot(-Umin,Amin,'r',-Umin,Amin-2*dAmin,'b:',-Umin,Amin+2*dAmin,'b:'), grid
    title('Minima - Choose a level where the estimate is stable'), 
    ylabel('Estimated mean exceedance, m'), xlabel('Lower threshold level, -u_{min}')
    
    Pout.Umin = Umin;
    Pout.Amin = Amin;
    Pout.Umax = Umax;
    Pout.Amax = Amax;
end

