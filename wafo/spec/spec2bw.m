function bw = spec2bw(S,fact)
%SPEC2BW Evaluates some spectral bandwidth and irregularity factors
%
% CALL:  bw = spec2bw(S,factors)
%
%        bw = vector of factors
%        S  = spectrum struct
%   factors = vector with integers, see below. (default [1])
%
% If input spectrum is of wave-number type, output are factors for
% corresponding 'k1D', else output are factors for 'freq'.
% Input vector 'factors' correspondence:
%    1 alpha=m2/sqrt(m0*m4)                        (irregularity factor)
%    2 eps2 = sqrt(m0*m2/m1^2-1)                   (narrowness factor)
%    3 eps4 = sqrt(1-m2^2/(m0*m4))=sqrt(1-alpha^2) (broadness factor)
%    4 Qp=(2/m0^2)int_0^inf f*S(f)^2 df            (peakedness factor)
% Order of output is the same as order in 'factors'
%
% Example:
%   S=demospec;
%   assert(spec2bw(S,[1 2 3 4]), [0.895608469201822, 0.230162951768626,...
%                                  0.444843196973911, 3.002022823753302], 1e-10)
% 

% References:
%

% Tested on: Matlab 5.3
% History: 
% Revised by jr 26.11.01
% - The variable vari was not assigned a value in the 
%   case of .type='freq'. Added an else statement in 
%   the second if sequence. 
% Revised by es 23.05.00
% - do not call spec2spec if already .type='freq'
% By es 23.09.1999

if nargin<2||isempty(fact)
  fact=1;
end

if isfield(S,'k')
  S=spec2spec(S,'k1d');
  vari='k';
elseif ~strcmpi(S.type,'freq')
  S=spec2spec(S,'freq');
  vari='w';
else
  vari = 'w';
end

m=spec2mom(S,4,[],0);
bw=zeros(size(fact));
for j=1:length(fact)
  switch fact(j)
    case 1
      bw(j)=m(3)/sqrt(m(1)*m(5));
    case 2
      bw(j)=sqrt(m(1)*m(3)/m(2)^2-1);
    case 3
      bw(j)=sqrt(1-m(3)^2/m(1)/m(5));
    case 4
      f = S.(vari);
      bw(j)=2/m(1)^2*simpson(f,f(:).*S.S(:).^2);
    otherwise
      error('Factor outside range (1,...,4)');
  end
end

    
  
