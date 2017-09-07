function h=lsplot(LS,cum,norm,beta,plottype)
%LSPLOT  Plot load spectra.
%
% CALL:  lsplot(LS)
%        h = lsplot(LS,cum,norm,beta,plottype)
%
% Input:
%   LS     = Load spectra.                            [N,2] OR {n}[N,2]
%   cum    = 0: Plot frequencies, 
%            1: Plot cumulative frequencies. (default)
%   norm   = 0: Use given frequencies, (default)
%            1: Normalize to relative frequencies.   
%   beta   = Plot damage spectrum with damage exponent beta. (default beta=0)
%   plottype = Plot load spectrum using commands: 
%              Cumulative plot: 'stairs' (default), or 'plot',
%              Frequency plot: 'bar' (default), 'stem', 'stairs', or 'plot'.
%
% Output:
%   h      = Handles to ploted lines.                        [n,1]
%
% The variable LS contains either one load spectrum stored in a matrix, or
% a cell-array of load spectra.  Each spectrum is a two collumn matrix with
% the amplitudes in the first column, and the frequencies in the second. 
%
% See also  cmat2amp, plotlc

% Copyright (c) 2003 by Pär Johannesson

% Tested  on Matlab  6.1, 6.5
%
% History:
% Created by TS (Thomas Svensson) 04-May-2002
% Updated by PJ 26-Jan-2003
% Updated by PJ 27-Jan-2003
%   Changed name from RFCPLOT to LSPLOT
%   Added plot of damage spectrum.
% Updated by PJ 03-Jun-2003
%   Added input argument 'plottype'
%   Changed name from PLOTLS to LSPLOT
% Updated by PJ 05-Jun-2003
%   Now correct cumulative plot.
%   Changed to 'stem' instead of 'plot'
% Updated by PJ 31-Oct-2003
%   Now cumulative plot starts from 1 (or minimal frequency).
%   Added 'bar' and 'stem' for frequency plot.

% Check input arguments
ni = nargin;
no = nargout;
%error(nargchk(1,5,ni));
narginchk(1,5)
if ni < 2
    cum=[];
end
if ni < 3
    norm=[];
end
if ni < 4
    beta=[];
end
if ni < 5
    plottype=[];
end

% Default settings
if isempty(cum)
    cum=1;
end
if isempty(norm)
    norm=0;
end
if isempty(beta)
    beta=0;
end
if isempty(plottype)
    if cum
        plottype='stairs';
    else
        plottype='bar';
    end
end

if isnumeric(LS), LS={LS}; end
n = length(LS);  % Number of spectra
h = zeros(n,1);  % Line handles

ih = ishold;
for k = 1:n
    LS{k}(:,2) = LS{k}(:,1).^beta.*LS{k}(:,2); % Damage spectrum if beta>0
    if beta > 0 
        if k==1, % Normalize damage to 1 for first spectrum
            norm_beta =sum(LS{k}(:,2)); % 
        end
        LS{k}(:,2) = LS{k}(:,2)/norm_beta;
    end
    if norm
        LS{k}(:,2) = LS{k}(:,2)/sum(LS{k}(:,2));
    end
    if cum
        I = (LS{k}(:,2)>0); LS1 = LS{k}(I,:); 
        if strcmp(plottype,'stairs')
            min_x = min(1,min(LS1(:,2)));
            x = [sum(LS1(:,2))-[0; cumsum(LS1(1:end-1,2))]; min_x; 0];
            y = [0; LS1(:,1); LS1(end,1)];
            
            %LS1=[[0; LS1(:,1)] [LS1(:,2); 0] ];
            %LS1 = LS1(size(LS1,1):-1:1,:);
            %h(k) = stairs(cumsum(LS1(:,2)),LS1(:,1));
            h(k) = stairs(flipud(x),flipud(y));
        elseif strcmp(plottype,'plot')
            LS1=[[0; LS1(:,1); LS1(end,1)] [0; LS1(:,2); 0] ];
            LS1 = LS1(size(LS1,1):-1:1,:);
            h(k) = plot(cumsum(LS1(:,2)),LS1(:,1));
        else
            error(['Invalid plottype: ' plottype])
        end
    else
        if strcmp(plottype,'bar')
            h(k) = bar(LS{k}(:,1),LS{k}(:,2),1);
        elseif strcmp(plottype,'stem')
            stem(LS{k}(:,1),LS{k}(:,2));
        elseif strcmp(plottype,'plot')
            h(k) = plot(LS{k}(:,1),LS{k}(:,2));
        elseif strcmp(plottype,'stairs')
            h(k) = stairs(LS{k}(:,1),LS{k}(:,2));
        else
            error(['Invalid plottype: ' plottype])
        end
    end
    hold on,
end

if ~ih, hold off, end

if cum
    set(gca,'xscale','log');
end


if beta == 0
    if cum == 0
        ylabel('Frequency of cycles');
        xlabel('Amplitude');
    else
        xlabel('Cumulative frequency of cycles');
        ylabel('Amplitude');
    end
else
    if cum == 0
        ylabel(['Damage per cycle, beta=' num2str(beta)]);
        xlabel('Amplitude');
    else
        xlabel(['Cumulative damage per cycle, beta=' num2str(beta)]);
        ylabel('Amplitude');
    end
end

if no==0, clear h; end
