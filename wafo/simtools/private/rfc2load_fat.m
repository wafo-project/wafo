function [X,res,comb,f]=rfc2load_fat(f,res,num_cc)
%RFC2LOAD_FAT  Recontructs a load process given the frequency matrix (and residual).
%
% CALL:  X = rfc2load_fat(f,num_cc,residual);
%
%  where
%
%        X        = the reconstructed load,
%        f        = the frequency matrix for the rainflow count,
%        residual = the residual (optional input argument),
%                   if left out the program will generate a
%                   stationary residual using the frequency 
%                   matrix.
%        num_cc   = the expected number of cycles,
%                   if num_cc=[] the program will continue
%                   until the residual is empty,

%  Copyright 1993, Mats Frendahl & Igor Rychlik,
%  Dept. of Math. Stat., University of Lund.

% Copyright (c) 2004 by Pär Johannesson

% Tested  on Matlab  6.5
%
% History:
% Adapted to WAFO by PJ (Pär Johannesson) 16-Feb-2004
%   The function 'rfc2load' originally from FAT (Fatigue Analysis Toolbox)
%   FAT is a predecessor of WAFO
% Changed by PJ 19-Feb-2004

%%%%
% Check input arguments

ni = nargin;
%no = nargout;
%error(nargchk(1,3,ni));
narginchk(1,3)
if ni<2, res=[]; end
if ni<3, num_cc=[]; end

% Check if freq. matrix is interger matrix.
if ~(round(f)==f)
    warning('WAFO:RFCLOAD_FAT','The frequency matrix is not interger matrix.   The matrix will be scaled.')
    f=floor(100*f);
end

N=length(f);

% If empty residual, then calculate a stat. residual
% using the frequency matrix.
if isempty(res)
    res0=fr2res(f);
    % Don't give a proper residual, modify it!
    % Correction by PJ
    [RFC,RFC1,res] = tp2rfc(res0(:));  % Ger proper residual
    resWAFO = N-res+1;              % Convert to WAFO-def
    RFMres = dtp2rfm(resWAFO,N,'CS'); % Rainflow matrix of residual
    RFMres = fliplr(RFMres)';       % Convert to FAT-def
    f = f-RFMres;                   % Remove cycles in res from rainflow matrix
end

% Find the largest amplitude. Look up how many there are,
% generate this number of large amplitudes, reset the 
% frequency matrix to 0 for this kind of amplitude.
largest_res=sort([min(min(res)) max(max(res))]);
index_i=largest_res(1); index_j=N-largest_res(2)+1;
n_large_res=f(index_i,index_j); f(index_i,index_j)=0;

xtra_cc=[];
for i=1:n_large_res, xtra_cc=[xtra_cc largest_res]; end

r=res(:)';

% Find the larges amplitude in the residual and put the 
% extra large amplitudes here.
undone=1;  
for i=1:length(r)-1
  if undone
    if ( sum(largest_res==[r(i) r(i+1)])==2 )
      r=[r(1:i-1) xtra_cc r(i:length(r))];
      undone=0;
    elseif ( sum(fliplr(largest_res)==[r(i) r(i+1)])==2 )
      r=[r(1:i-1) fliplr(xtra_cc) r(i:length(r))];
      undone=0;
    end
  end
end


% If the number of requested cycles are not empty then find
% the number of extected runs to give this number of cycles.
if ~isempty(num_cc)
    % Norm frequency matrix to probability matrix.
    f=f/sum(sum(f));
        
    % Calculate the number of extected runs to get the requested number 
    % of cycles.
    num_cc_limit=2*num_cc+4;
    
    % Scale the frequency matrix to whole numbers and make it large to 
    % be stationary.
    f=floor(max([1000*N 1e6])*f);
end

% Calculated the combination matrix
comb=fr2comb(f,r);

% Initiate the residual.
X=r(1); 

again=1; num_runs=0;

while again
    
    num_runs=num_runs+1;
        
    h=0;
    
    if r(1)>r(2)
        % disp('upslope')
        
        if r(1)>(r(2)+1)
            m=r(1); M=r(2); j=N-m+1; i=m-1:-1:M+1; 
            combinations=comb(i,j);
            freq=f(i,j);
            prob=comb2pro(combinations,freq);
            h=min(find(rand<cumsum(prob)))-1;
        end
        
        if h>0
            C=m-h;
            %i=C:m; 
            if f(C,j)>0
                f(C,j)=f(C,j)-1;
                i=(C+1):r(1);
                j=N+1-r(1);
                comb(i,j)=comb(i,j)-1;
                r=[m-1 C r]; 
            else
                warning('Warning 1')
                break
                %r=[r(1)-1 r(2:length(r))];
            end
        else
            j=N+1-r(1);
            if r(1)>(r(2)+1)
                i=(r(2)+1):(r(1)-1);
                comb(i,j)=comb(i,j)-1;
            end
            r(1)=r(1)-1;
            
        end
        
    else
        
        % disp('downslope')
        if r(1)<r(2)-1
            m=r(2); M=r(1); i=M; j=M+1:m-1; j=N-j+1;
            combinations=comb(i,j)';
            freq=f(i,j)';
            prob=comb2pro(combinations,freq);
            h=min(find(rand<cumsum(prob)))-1;
        end   
        
        if h>0
            j=M:M+h; j=N-j+1; 
            c0=M+h;
            c=N+1-c0;
            if f(i,c)>0
                f(i,c)=f(i,c)-1;
                i=r(1);
                j=N+1-((i+h-1:-1:r(1)));
                comb(i,j)=comb(i,j)-1;
                r=[M+1 c0 r];
            else
                warning('Warning 2'), break
                %r=[r(1)+1 r(2:length(r))];
            end
        else
            i=r(1);
            j=N+1-(r(1):(r(2)-1));
            comb(i,j)=comb(i,j)-1;
            r(1)=r(1)+1;
                                    
        end
    end
    
    if r(1)==r(2), r=r(2:length(r)); end  
    
    if (length(r)>=2)
        again=1;
    else
        disp('   The residual is empty. Program will terminate.')
        again=0;
    end
    
    X=[X r(1)];
    
    n_tmp=length(X);
    
    % This is an equivalent of shave, but faster
    
    if n_tmp>2
        if (X(n_tmp-2)<X(n_tmp-1))&&(X(n_tmp-1)<X(n_tmp))
            X=[X(1:n_tmp-2) X(n_tmp)];
        elseif (X(n_tmp-2)>X(n_tmp-1))&&(X(n_tmp-1)>X(n_tmp))
            X=[X(1:n_tmp-2) X(n_tmp)];
        end
    end
    
    comb=fliplr(triu(fliplr(comb),1));
    
    % If the number of requested cycles are not empty then check
    % if the number of step (length of load) exceeds the expected 
    % number.
    if ~isempty(num_cc)
        if length(X)>num_cc_limit
            disp('   The number of requested cycles has been reached.')
            disp('   Program will terminate.')
            again=0;
        end
    end
    
end

X = X';
