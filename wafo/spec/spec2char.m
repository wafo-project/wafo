function [ch,R1,chtext,R]=spec2char(S,fact,T)
%SPEC2CHAR  Evaluates spectral characteristics and their covariance
%
% CALL: [ch R chtext] = spec2char(S,fact,T)
%
%       ch = vector of spectral characteristics
%       R  = matrix of the corresponding covariances given T
%   chtext = a cellvector of strings describing the elements of ch, see example.
%       S  = spectral struct with angular frequency
%     fact = vector with factor integers or a string or
%            a cellarray of strings, see below.(default [1])
%       T  = recording time (sec) (default 1200 sec = 20 min)
%
% If input spectrum is of wave number type, output are factors for
% corresponding 'k1D', else output are factors for 'freq'.
% Input vector 'factors' correspondence:
%    1 Hm0   = 4*sqrt(m0)                              Significant wave height
%    2 Tm01  = 2*pi*m0/m1                              Mean wave period
%    3 Tm02  = 2*pi*sqrt(m0/m2)                        Mean zero-crossing period
%    4 Tm24  = 2*pi*sqrt(m2/m4)                        Mean period between maxima
%    5 Tm_10 = 2*pi*m_1/m0                             Energy period
%    6 Tp    = 2*pi/{w | max(S(w))}                    Peak period  
%    7 Ss    = 2*pi*Hm0/(g*Tm02^2)                     Significant wave steepness
%    8 Sp    = 2*pi*Hm0/(g*Tp^2)                       Average wave steepness
%    9 Ka    = abs(int S(w)*exp(i*w*Tm02) dw ) /m0     Groupiness parameter
%   10 Rs    = (S(0.092)+S(0.12)+S(0.15)/(3*max(S(w))) Quality control parameter
%   11 Tp1   = 2*pi*int S(w)^4 dw                      Peak Period (robust estimate for Tp) 
%              ------------------
%              int w*S(w)^4 dw 
%
%   12 alpha = m2/sqrt(m0*m4)                          Irregularity factor
%   13 eps2  = sqrt(m0*m2/m1^2-1)                      Narrowness factor
%   14 eps4  = sqrt(1-m2^2/(m0*m4))=sqrt(1-alpha^2)    Broadness factor
%   15 Qp    = (2/m0^2)int_0^inf w*S(w)^2 dw           Peakedness factor
%
% Order of output is same as order in 'factors'
% The covariances are computed with a Taylor expansion technique
% and is currently only available for factors 1, 2, and 3. Variances
% are also available for factors 4,5,7,12,13,14 and 15 
% 
% Quality control:
%  Critical value for quality control parameter Rs is Rscrit = 0.02
%  for surface displacement records and Rscrit=0.0001 for records of
%  surface acceleration or slope. If Rs > Rscrit then probably there 
%  are something wrong with the lower frequency part of S.
%
%  Ss may be used as an indicator of major malfunction, by checking that
%  it is in the range of 1/20 to 1/16 which is the usual range for
%  locally generated wind seas. 
%
% Examples:
%   S      = demospec;
%   [ch R,txt] = spec2char(S,[1 2 3]);    % fact a vector of integers
%   cat(2,char(txt),repmat(' = ',3,1),num2str(ch')) 
%   ch = num2cell(ch);
%   ch0 = cell2struct(ch,txt,2);          % Make a structure 
%   [txt{2,:}]=deal(','); txt{2,end}=' ';
%   eval(['[' txt{:} '] = deal(ch{:})'])  % Assign values to variables
%   ch1    = spec2char(S,'Ss');           % fact a string
%   ch2    = spec2char(S,{'Tp','Tp1'});   % fact a cellarray of strings
%
% See also  dspec2char, spec2bw, spec2mom, simpson

% References:
% Krogstad, H.E., Wolf, J., Thompson, S.P., and Wyatt, L.R. (1999)
% 'Methods for intercomparison of wave measurements'
% Coastal Enginering, Vol. 37, pp. 235--257
%
% Krogstad, H.E. (1982)
% 'On the covariance of the periodogram'
% Journal of time series analysis, Vol. 3, No. 3, pp. 195--207
%
% Tucker, M.J. (1993)
% 'Recommended standard for wave data sampling and near-real-time processing'
% Ocean Engineering, Vol.20, No.5, pp. 459--474
%
% Young, I.R. (1999)
% "Wind generated ocean waves"
% Elsevier Ocean Engineering Book Series, Vol. 2, pp 239

% Tested on: Matlab 7
% History: 
% revised pab
% updated to matlab 7 syntax
% revised pab jan2004
%  -added todo comments
% revised pab 25.06.2001
% -added chtext to output+more examples
% revised pab 07.07.2000
%  -fixed a bug for variance of Tm02
%  - added variance for Tm24
% revised pab 22.06.2000
%  - added alpha, eps2,eps4,Qp
%  - added the possibility that fact is a string or a cellarray of strings 
% revised pab 23.05.2000
%  - added ttspec call
%  revised by es 23.05.00, do not call spec2spec if already .type='freq'
% Revised pab 06.03.2000
% -updated header info
%by pab 16.02.2000
  
  
% TODO % Need more checking on computing the variances for Tm24,alpha, eps2 and eps4 
% TODO % Covariances between Tm24,alpha, eps2 and eps4 variables are also needed

tfact = {'Hm0', 'Tm01', 'Tm02', 'Tm24', 'Tm_10','Tp','Ss', 'Sp', 'Ka', ...
      'Rs', 'Tp1','alpha','eps2','eps4','Qp'} ;

if nargin<2||isempty(fact)
  nfact = 1;
elseif iscell(fact)||ischar(fact)
  if ischar(fact), fact = {fact}; end
  N     = length(fact(:));
  nfact = zeros(1,N);
  ltfact = char(lower(tfact));
  for ix=1:N,
    ind = strmatch(lower(char(fact(ix))),ltfact,'exact');
    if length(ind)==1,
      nfact(ix)=ind;
    else
      error(['Not a valid factor: ' fact{ix}])
    end
  end 
else
  nfact = fact;
end
if any(nfact>15 | nfact<1)
  error('Factor outside range (1,...,15)')
end

if nargin<3||isempty(T)
  T = 1200; % recording time default 1200 secs (=20 minutes)
end

if isfield(S,'k')
  S=spec2spec(S,'k1d');
  vari='k';
else 
  if ~strcmpi(S.type,'freq'),
    S = spec2spec(S,'freq');
  end
  S = ttspec(S,'w');
  vari = 'w';
end

f  = S.(vari);
f  = f(:);
S1 = S.S(:);
%m=spec2mom(S,4,[],0)./[ (2*pi).^[0:4] ]; % moments corresponding to freq
% in Hz
m = simpson(f,[S1 f.*S1 f.^2.*S1 f.^3.*S1 f.^4.*S1])./ ((2*pi).^(0:4) );

ind  = find(f>0);
m(6) = simpson(f(ind),S1(ind)./f(ind))*2*pi;  % = m_1
m_10 = simpson(f(ind),S1(ind).^2./f(ind))*(2*pi)^2/T;    % = COV(m_1,m0|T=t0)
m_11 = simpson(f(ind),S1(ind).^2./f(ind).^2)*(2*pi)^3/T; % = COV(m_1,m_1|T=t0)


%      Hm0        Tm01        Tm02             Tm24         Tm_10
Hm0  = 4*sqrt(m(1)); 
Tm01 = m(1)/m(2); 
Tm02 = sqrt(m(1)/m(3)); 
Tm24 = sqrt(m(3)/m(5)); 
Tm_10= m(6)/m(1);

Tm12 = m(2)/m(3);

g    = gravity;
[maxS ind] = max(S1);
Tp   = 2*pi/f(ind);                                   % peak period /length
Ss   = 2*pi*Hm0/g/Tm02^2;                             % Significant wave steepness
Sp   = 2*pi*Hm0/g/Tp^2;                               % Average wave steepness 
Ka   = abs(simpson(f,S1.*exp(sqrt(-1)*f*Tm02)))/m(1); % groupiness factor

% Quality control parameter 
% critical value is approximately 0.02 for surface displacement records
% If Rs>0.02 then there are something wrong with the lower frequency part 
% of S.
Rs   = sum(interp1(f,S1,[0.0146 0.0195 0.0244]*2*pi,'linear'))/3/maxS; 
Tp2  = 2*pi*simpson(f,S1.^4)/simpson(f,f.*S1.^4);


alpha1 = Tm24/Tm02;                 % m(3)/sqrt(m(1)*m(5));
eps2   = sqrt(Tm01/Tm12-1);         % sqrt(m(1)*m(3)/m(2)^2-1);
eps4   = sqrt(1-alpha1^2);          % sqrt(1-m(3)^2/m(1)/m(5));
Qp     = 2/m(1)^2*simpson(f,f.*S1.^2);




ch = [Hm0 Tm01 Tm02 Tm24 Tm_10 Tp Ss Sp Ka Rs Tp2 alpha1 eps2 eps4 Qp];


% Select the appropriate values
ch     = ch(nfact);
chtext = tfact(nfact);

if nargout>1,
  % covariance between the moments:
  %COV(mi,mj |T=t0) = int f^(i+j)*S(f)^2 df/T
  ONE = ones(size(f));
  S2 = S1.^2;
  mij = simpson(f,[ONE f f.^2  f.^3 f.^4 f.^5 f.^6 f.^7 f.^8].*S2(:,ones(1,9)))/T./( (2*pi).^(-1:7) );
%   mij = simpson(f,[S1.^2 f.*S1.^2 (f.*S1).^2  f.^3.*S1.^2 (f.^2.*S1).^2 ...
%	f.*(f.^2.*S1).^2 (f.^3.*S1).^2 f.*(f.^3.*S1).^2 (f.^4.*S1).^2])/T./[ (2*pi).^[-1:7] ]

% and the corresponding variances for
%{'hm0', 'tm01', 'tm02', 'tm24', 'tm_10','tp','ss', 'sp', 'ka', 'rs', 'tp1','alpha','eps2','eps4','qp'}  
  R = [4*mij(1)/m(1) ...
	mij(1)/m(2)^2-2*m(1)*mij(2)/m(2)^3+m(1)^2*mij(3)/m(2)^4 ...
	0.25*(mij(1)/(m(1)*m(3))-2*mij(3)/m(3)^2+m(1)*mij(5)/m(3)^3) ...
	0.25*(mij(5)/(m(3)*m(5))-2*mij(7)/m(5)^2+m(3)*mij(9)/m(5)^3) ...
	m_11/m(1)^2+(m(6)/m(1)^2)^2*mij(1)-2*m(6)/m(1)^3*m_10,...
	NaN,...
	(8*pi/g)^2*(m(3)^2/(4*m(1)^3)*mij(1)+mij(5)/m(1)-m(3)/m(1)^2*mij(3)),...
	NaN*ones(1,4),...
	m(3)^2*mij(1)/(4*m(1)^3*m(5))+mij(5)/(m(1)*m(5))+mij(9)*m(3)^2/(4*m(1)*m(5)^3)-...
	m(3)*mij(3)/(m(1)^2*m(5))+m(3)^2*mij(5)/(2*m(1)^2*m(5)^2)-m(3)*mij(7)/m(1)/m(5)^2,...
	(m(3)^2*mij(1)/4+(m(1)*m(3)/m(2))^2*mij(3)+m(1)^2*mij(5)/4-m(3)^2*m(1)*mij(2)/m(2)+...
        m(1)*m(3)*mij(3)/2-m(1)^2*m(3)/m(2)*mij(4))/eps2^2/m(2)^4,...
	(m(3)^2*mij(1)/(4*m(1)^2)+mij(5)+m(3)^2*mij(9)/(4*m(5)^2)-m(3)*mij(3)/m(1)+....
	m(3)^2*mij(5)/(2*m(1)*m(5))-m(3)*mij(7)/m(5))*m(3)^2/(m(1)*m(5)*eps4)^2,...
	NaN];
 
  % and covariances by a taylor expansion technique:
  % Cov(Hm0,Tm01) Cov(Hm0,Tm02) Cov(Tm01,Tm02)
  S0 = [ 2/(sqrt(m(1))*m(2))*(mij(1)-m(1)*mij(2)/m(2)),...
	1/sqrt(m(3))*(mij(1)/m(1)-mij(3)/m(3)),......
	1/(2*m(2))*sqrt(m(1)/m(3))*(mij(1)/m(1)-mij(3)/m(3)-mij(2)/m(2)+m(1)*mij(4)/(m(2)*m(3)))];
  tmp = NaN;
  R1  = tmp(ones(15,15));
  
  for ix=1:length(R), 
    R1(ix,ix) = R(ix);  
  end
  
  
  R1(1,2:3)   = S0(1:2);
  R1(2,3)     = S0(3);
  for ix = 1:2, %make lower triangular equal to upper triangular part
    R1(ix+1:3,ix) = R1(ix,ix+1:3).';
  end
  
  R = R(nfact);
  R1= R1(nfact,nfact);
end

 % Needs further checking:
 % Var(Tm24)= 0.25*(mij(5)/(m(3)*m(5))-2*mij(7)/m(5)^2+m(3)*mij(9)/m(5)^3) ...



