function [b]= subsref(self,index)
% SUBSREF Define field name indexing or methods to apply on the objects 
%
%  CALL b = subsref(self,index)
%
% Examples
%  S = mkjonswap();
%  S.Hm0
%  ft = S.freqtype;
%  plot(S.(ft),S.S)            % Plot spectrum
%  S.plot                      % Alternative call
%  argMax = S.argMax;          % = S.wdata.argMax
%  argMax = S.wdata.argMax;
%
%  S2 = S; S2(2) = specdata(torsethaugen);
%  S2.plot       % Plot both spectra in separate figures
%
%  am = {S2.argMax}            % cellarray of argMax values
%  am1 = {S2(1:2).argMax}      % NB!: Matlabbug returns only the first value
%  [am2{1:2}] = S2(1:2).argMax % Use this call instead!
%  
% See also specdata/get


% TODO % This is slow. Optimize this for speed.
% TODO % am1 = {S2(1:2).argMax} is a Matlabbug returns only the first value

%History
%By pab 2016

  N = length(index);
  subs1 = index(1).subs;
  switch index(1).type
    case '()'
      wn = subs1{1};
      b = jonswap(self, wn);
    case '{}'
      error('WAFO:MKJONSWAP','??? Cell contents reference from a non-cell array object.')
    case '.'
      b = self.(subs1);
    otherwise
        error('WAFO:MKJONSWAP','???? Reference to non-existent field or method %s', subs1)
  end % switch
end


%% sub function
  function S = jonswap(options, w)
  %JONSWAP spectral density
  % 
  % CALL      S = jonswap(w)
  %     options = jonswap('options')
  %
  

    if strncmpi(w,'options',1)
      S = struct(options);
    else
      if (options.Hm0>0)
      
        N   = options.N;
        M   = options.M;
        wp  = 2*pi/options.Tp;
        wn  = w/wp;
        Ag  = options.Ag;
        Hm0 = options.Hm0;
        Gf  = peakEnhancementFactor(wn,options);
        S   = ((Hm0/4)^2/wp.*Ag).*Gf.*ggamspec(wn,N,M);
      else
        S = zeros(size(w));
      end
    end
  end % jonswap


