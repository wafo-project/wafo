function ind=findrfc(tp,h)
%FINDRFC Finds indices to rainflow cycles of a sequence of TP.
%
% CALL:  RFC_ind = findrfc(TP,h);
%
%        TP  = vector of turningpoints (NB! Only values, not sampled times)
%
%    RFC_ind = indices to the rainflow cycles of the original sequence TP.
%
%         h  = a threshold, must be larger than zero. 
%              if h>0, then all rainflow cycles with height 
%              smaller than  h  are removed.
%
%  This function is not implemented as a matlab function; instead, a 
%  mex file (originally written in C) is utilized.
% 
%  Example:
%   x = load('sea.dat'); 
%   tp = dat2tp(x); 
%   ind = findrfc(tp(:,2),0.3); 
%   waveplot(x,tp(ind,:),1,1);
%
%   assert(length(ind), 1008);
%   assert(ind(1:5)', [1, 2, 7, 8, 9]);
%   close all;
% 
%  See also  tp2rfc, dat2tp,  rfcfilter. 
  
% This is a modified version of rfcfilt (found in WAT), which is about 20 
% to 30 times faster than rfcfilt (on a PentiumII 233 MHz  
% with 32 MB ram and Matlab 5.0 under Linux). The reason is
% that this version does not save TP to disk. Instead it passes 
% the arguments directly to the executeable file. 
% However, this solution requires different input and returns
% indices to the rfc instead of the rfc itself.
% It also ignores the first turningpoint if that is a maximum and 
% starts on the first minimum when finding the sequence of rfc. 

% Tested on Matlab 6.1,6.0, 5.2
% History:
% revised pab may 2007
% - translated c-code to matlab.
% revised pab 10.08.2003
% - fixed a bug in the example  
% Revised by pab 24.07.1999


%disp('FINDRFC is not implemented as a m-function')
%disp('                   compile the mexfile findrfc.c before you try again.')
%error('findrfc error')

%error(nargchk(2,2,nargin))
narginchk(2,2)   
ind = [];
isFirstAmax = tp(1)>tp(2);
Tstart=1;
if isFirstAmax % if first is a max*/
  tp(1) = []; % ignore the first max*/
  Tstart=2;
end

n = length(tp);
NC=floor(n/2);

 
if (NC<1)
  return % No RFC cycles*/
end
   
dtp = diff(tp);
if any(dtp(1:n-2).*dtp(2:n-1)>=0)
  warning('WAFO:findrfc','This is not a sequence of turningpoints! Exiting.')
  return
end
   

ix = 0;
   
ind = zeros(size(tp));

%for ii=0:NC-1
for ii=0:NC-2
  
  Tmi = Tstart+2*ii;
  Tpl = Tstart+2*ii+2;
  xminus = tp(2*ii+1);
  xplus = tp(2*ii+2+1);

  if(ii~=0)
    j=ii-1;
    while((j>=0) && (tp(2*j+1+1)<=tp(2*ii+1+1)))
      if ( tp(2*j+1)<xminus) 
        xminus=tp(2*j+1);
        Tmi=Tstart+2*j;
      end %} /*if */
      j = j-1; %j--;
    end %/*while j*/
  end % } /*if ii */
  
  if ( xminus >= xplus)
    if ( (tp(2*ii+1+1)-xminus) >= h)
      ix = ix+1;
      ind(ix)=Tmi;
      ix = ix+1;% ix++;
      ind(ix) = (Tstart+2*ii+1);
	   
    end %} /*if*/
    %goto L180;
    %end %if  }
  else
      
    j=ii+1;
    gotoL170 = false;
    %while(j<NC )
    while(j<NC-1 )
       gotoL170 = (tp(2*j+1+1) >= tp(2*ii+1+1));
      if gotoL170 
        break
      end
      if( tp(2*j+2+1) <= xplus)
        xplus=tp(2*j+2+1);
        Tpl = (Tstart+2*j+2);
      end %if/
      j = j +1; %++;
    end %} /*while*/
    
    if gotoL170
      % L170:
      if (xplus <= xminus )
        if ( (tp(2*ii+1+1)-xminus) >= h)
          ix = ix+1;
          ind(ix)=Tmi;
          ix = ix+1;
          ind(ix) = (Tstart+2*ii+1);
        end %} /*if*/
        %/*goto L180;*/
        %   }
      elseif ( (tp(2*ii+1+1)-xplus) >= h)
        ix = ix+1;
        ind(ix)=(Tstart+2*ii+1);
        ix = ix+1;
        ind(ix)=Tpl;
      end %} /*elseif*/      
    elseif ((tp(2*ii+1+1)-xminus) >= h)
      ix = ix+1;
      ind(ix)=Tmi;
      ix = ix +1;
      ind(ix)=(Tstart+2*ii+1);
    end %} /*if*/
    %goto L180;
  end % %%%%
  % L180:
  % iy=ii;
end %}  /* for ii */
%ind(ix+1:end) = [];

ind = sort(ind(1:ix));
%return ix;

