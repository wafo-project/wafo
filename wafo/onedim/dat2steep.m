function [S , H,AC1,AT1,TFRONT1,TREAR1,z_ind2,xn2]=dat2steep(xx,rate,method)
%DAT2STEEP Extracts waveheights and steepnesses from data.
%
% CALL:  [S, H,Ac,At,Tcf,Tcb, z_ind,yn] = dat2steep(xn,rate,method);
%
%   S, H = Steepness and the corresponding wave height according to method
%  Ac,At = crest and trough amplitude, respectively
%   Tcf,
%    Tcb = Crest front and crest (rear) back period, respectively
%  z_ind = indices to the zero-crossings (d,u) of the defining  
%          trough to trough waves (tw). If M>1 then 
%          z_ind=[N1 z_ind1 N2 z_ind2 ...NM z_indM] where 
%          Ni = length(z_indi) and z_indi are the indices to 
%          the zero-crossings of xi, i=1,2...M.
%
%     yn = interpolated signal
%
%     xn = [ti x1 x2 ... xM], where 
%          ti = time and x1 x2 ... xM are M column vectors of 
%          sampled surface elevation.
%
%   rate = 1,2,3..., interpolation rate  
%          no interpolation done before extracting the
%          parameters if less than one. Interpolates
%          with spline if greater than one.
%                          
% method = 0 max(Vcf, Vcb) and corresponding wave height Hd or Hu in H
%          1 crest front (rise) speed (Vcf) in S and wave height Hd in H. (default)
%         -1 crest back (fall) speed (Vcb) in S and waveheight Hu in H.
%          2 crest front steepness in S and the wave height Hd in H.
%         -2 crest back steepness in S and the wave height Hu in H.
%          3 total wave steepness in S and the wave height Hd in H
%            for zero-downcrossing waves.
%         -3 total wave steepness in S and the wave height Hu in H.
%            for zero-upcrossing waves.
%
% The parameters are calculated as follows:
%  Crest front speed (velocity) = Vcf = Ac/Tcf
%  Crest back speed  (velocity) = Vcb = Ac/Tcb
%  Crest front steepness  =  2*pi*Ac./Td/Tcf/g
%  Crest back steepness   =  2*pi*Ac./Tu/Tcb/g
%  Total wave steepness (zero-downcrossing wave) =  2*pi*Hd./Td.^2/g
%  Total wave steepness (zero-upcrossing wave)   =  2*pi*Hu./Tu.^2/g
%   
% The definition of g, Ac,At, Tcf, etc. are given in gravity, wavedef,
% ampdef, and perioddef. 
%  
% Example:
%  dt = 0.4;
%  xs = spec2sdat(specinterp(jonswap,dt),6000); rate=8; method=1;
%  [S,H] = dat2steep(xs,rate,method);
%  plot(S,H,'.'),xlabel('Vcf [m/s]'),ylabel('Hd [m]')
%
% See also  wavedef, ampdef, perioddef, interp1, dat2tc

% Tested on:Matlab 5.3, 5.2, 5.1
% History:
% revised pab Feb2004  
% revised pab 01.12.2002
% -removed disp statement and replaced with call to waitbar.
% revised pab 13.06.2001
% - changed method +/-4 to +/-3
% - added call to ecross => improved accuracy in the zero-crossing
%   period calculations
% revised pab 03.04.2001
% - changed order of methods (hopefully more logical)
% - added example
% revised pab 28.11.2000
% -fixed a bug when rate=1 and M>2
% revised pab 07.03.2000
%  - added method 4
% revised pab 08.02.2000
%  - added Ac,At,Tcf, Tcb to output
% by pab 11.11.98

%    
  
  error(nargchk(1,3,nargin))
if nargin<3||isempty(method),  
  method=1; % want crestfrontvelocity
end

[S,H,z_ind2,AC1,AT1,TFRONT1,TREAR1]=deal([]); % Initialize to []
dT    = xx(2,1)-xx(1,1);
[N M] = size(xx);
g     = gravity;  % acceleration of gravity
interpolate = 0;

if nargin<2||isempty(rate)||(rate<=1),  % no interpolation
  xn=xx(:,1:2); 
else % interpolate with spline
  dT  = dT/rate;
  ti  = (xx(1,1):dT:xx(N,1))';
  interpolate=1;
  xn  = zeros(length(ti),2);
  xn(:,1) = ti;
  if nargout>=7
    xn2 = zeros(length(ti),M);
    xn2(:,1) = ti;
  end
end

%h9 = waitbar(0,'Please wait... for dat2steep to finish.');
for ix=2:M
  waitbar((ix-1)/(M-1)); %,h9)
  if interpolate,
    %disp(['   ...interpolating column ' int2str(ix)])
    %    xn=[ti ;  interp1(xx(:,1),xx(:,2),ti,'*linear') ]'; 
    xn(:,2)=interp1(xx(:,1),xx(:,ix),ti,'*spline'); 
    %    xn=[interp(xx(:,1),rate)';interp(xx(:,2),rate)' ]';
    %    disp(' finished')
     if nargout>7
       xn2(:,ix)=xn(:,2);
     end
   else
     xn(:,2)=xx(:,ix);
  end
 % disp('   ...extracting parameters')
  [tc tc_ind,z_ind]=dat2tc(xn,0,'tw'); % finding trough to trough waves
  if nargout>6
    if M==2,
      z_ind2=z_ind; % indices to zerocrossings of xn
    else
      z_ind2 = [z_ind2; length(z_ind); z_ind];
    end
  end
  % crest amplitude
  AC=tc(2:2:end,2);
  % trough  amplitude
  AT=-tc(1:2:end,2);

  if (0<= method && method <=2)||nargout>4,
    % time between zero-upcrossing and  crest  [s]
    tu     = ecross(xn(:,1),xn(:,2), z_ind(2:2:(end-1)),0);
    TFRONT = tc(2:2:end,1)-tu;
    %TFRONT=tc(2:2:end,1)-xn(z_ind(2:2:(end-1)),1);
    TFRONT((TFRONT==0))=dT; % avoiding division by zero
  end
  if (0 >= method && method>=-2)||nargout>5,
    % time between  crest and zero-downcrossing [s]
    td    = ecross(xn(:,1),xn(:,2), z_ind(3:2:end),0);
    TREAR = td-tc(2:2:end,1);
    %TREAR=xn(z_ind(3:2:end),1)-tc(2:2:end,1);
    TREAR((TREAR==0))=dT;  % avoiding division by zero
  end   
 
  switch method
    case 0,   % max(Vcf, Vcr) and the corresponding wave height Hd or Hu in H
      HU      = AC+AT(2:end);
      [T ind] = min([TFRONT TREAR],[],2);
      
      H2      = AC+AT(1:end-1);
      S       = [S ;AC./T];
      ind     = (ind==2);
      H2(ind) = HU(ind);
      H       = [H;H2];
    case 1,  % extracting crest front velocity [m/s] and  
             % Zero-downcrossing wave height [m]
      H = [H;AC+AT(1:end-1)] ; % Hd
      S = [S ; AC./TFRONT];
      
    case -1,% extracting crest rear velocity [m/s] and  
           % Zero-upcrossing wave height [m]
      H = [H ; AC+AT(2:end)]; 
      S = [S ; AC./TREAR];
      
    case 2,  % crest front steepness in S and the wave height Hd in H.
      H = [H;AC+AT(1:end-1) ];
      T = diff(ecross(xn(:,1),xn(:,2), z_ind(1:2:end),0));
      %T =xn(z_ind(3:2:end),1)-xn(z_ind(1:2:(end-2)),1);
      S = [S ; 2*pi*AC./T./TFRONT/g];
    case -2,  % crest back steepness in S and the wave height Hu in H.
      H = [H;AC+AT(2:end) ];
      T = diff(ecross(xn(:,1),xn(:,2), z_ind(2:2:end),0));
      %T=xn(z_ind(4:2:end),1)-xn(z_ind(2:2:(end-1)),1);
      S = [S ; 2*pi*AC./T./TREAR/g];
    case 3,   % total steepness in S and the wave height Hd in H
              % for zero-doewncrossing waves.
      H = [H;AC+AT(1:end-1) ];
      T = diff(ecross(xn(:,1),xn(:,2), z_ind(1:2:end),0));
      %T=xn(z_ind(3:2:end),1)-xn(z_ind(1:2:(end-2)),1); % Period zero-downcrossing waves
      S = [S ; 2*pi*H./T.^2/g];
    case -3,   % total steepness in S and the wave height Hu in H for
               % zero-upcrossing waves.
      H = [H;AC+AT(2:end) ];
      T = diff(ecross(xn(:,1),xn(:,2), z_ind(2:2:end),0));
      %T=xn(z_ind(4:2:end),1)-xn(z_ind(2:2:(end-1)),1); % Period zero-upcrossing waves
      S = [S ; 2*pi*H./T.^2/g]; 
    otherwise,  error('unknown option')
   end
   if nargout>2
     AC1=[AC1;AC];
     if nargout>3
       AT1=[AT1;AT];
       if nargout>4
	 TFRONT1=[TFRONT1;TFRONT];
	 if nargout>5
	   TREAR1=[TREAR1;TREAR];
	 end
       end
     end
   end
    
   if 0,
     ind=(AC<0.5);
     V(ind)=[];
     H(ind)=[];
   end
end 
%close(h9)
return

