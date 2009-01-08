function [m,mtext] = spec2mom(S,nr,vari,even)
%SPEC2MOM Calculates spectral moments from spectrum
%
% CALL:  [m,mtext] = spec2mom(S,nr,variables,even)
%
%   m    = vector of moments
%   mtext= a cell array of strings describing the elements of m, see below
%   S    = spectrum (struct)
%   nr   = order of moments (default 2, maximum 4)
%   variables = variables in model, optional when two-dim.spectrum,
%               string with x and/or y and/or t  (default 'xt')
%   even = 0 for all moments, 1 for only even orders (default 1)
%
% Calculates spectral moments of up to order four by use of
% Simpson-integration.
%
%       //
% m_jkl=|| k1^j*k2^k*w^l S(w,th) dw dth 
%       //
%
% where k1=w^2/gravity*cos(th-phi),  k2=w^2/gravity*sin(th-phi)
% and phi is the angle of the rotation in S.phi. If the spectrum 
% has field .g, gravity is replaced by S.g.
%  
% The strings in output mtext have the same position in the cell array
% as the corresponding numerical value has in output m
% Notation in mtext: 'm0' is the variance,
%                    'mx' is the first-order moment in x,
%                   'mxx' is the second-order moment in x,
%                   'mxt' is the second-order cross moment between x and t,
%                 'myyyy' is the fourth-order moment in y
%                         etc.
% For the calculation of moments see Baxevani et al.
% Example:
%      S=demospec('dir')
%      [m,mtext]=spec2mom(S,2,'xyt')

% References
%  Baxevani A. et al. (2001) 
%  Velocities for Random Surfaces 
%
%
% Tested on: Matlab 6.0
% Tested on: Matlab 5.3
% History:
% Revised by I.R. 04.04.2001: Introducing the rotation angle phi.
% Revised by A.B. 23.05.2001: Correcting 'mxxyy' and introducing
% 'mxxyt','mxyyt' and 'mxytt'.
% Revised by A.B. 21.10.2001: Correcting 'mxxyt'.
% Revised by A.B. 21.10.2001: Adding odd-order moments.
% By es 27.08.1999


if nargin<2||isempty(nr)
  nr=2;
end
if nargin<4||isempty(even)
  even=1;
end
even=logical(even);

if strcmpi(S.type,'freq')|| strcmpi(S.type,'enc') ||...
      strcmpi(S.type,'k1d') %i.e. one-dim spectrum
  if isfield(S,'w')
    f=S.w(:);
    vari='t';
  elseif isfield(S,'k')
    f=S.k(:);
    vari='x';
  else
    f=2*pi*S.f(:);
    S.S=S.S/2/pi;
    vari='t';
  end
  S1=abs(S.S(:));
  m=simpson(f,S1);
  mtext={'m0'};
  if nr>0
    if ~even
      m=[m simpson(f, f.*S1)];
      mtext(end+1)={'mi'};
    end
    if nr>1
      m=[m simpson(f,f.^2.*S1)];
      mtext(end+1)={'mii'};
      if nr>2
        if ~even
          m=[m simpson(f,f.^3.*S1)];
          mtext(end+1)={'miii'};
        end
        if nr>3
          m=[m simpson(f,f.^4.*S1)];
          mtext(end+1)={'miiii'};
        end
      end
    end
    mtext(end+1)={strcat(' where i=',vari)};
  end






else %two-dim spectrum
  if (nargin<3 || isempty(vari)) && nr<=1;
    vari='x';% set default value
    Nv=1;
    elseif nargin<3||isempty(vari);
    vari='xt';
    Nv=2;
    else% secure the mutual order ('xyt')
    vari=sort(lower(vari));
    Nv=length(vari);
    
    if ( (vari(1)=='t') && (Nv>1)),
      vari=[vari(2:Nv), 't'];
    end 
  end
  S1=S;
  %%%%to be removed
  %%%%%%%
  %if Nv==1
   % switch vari
    %  case 'x'
        %S1=spec2spec(S,'k1d');
     % case 'y'
%       S1=spec2spec(S,'k1d',-pi/2); % funkar detta??? -pi/2 ???
    %  case 't'
%       S1=spec2spec(S,'freq');

%end
 %   [m,mtext]=spec2mom(S1,nr,vari,even);
 % else % Nv>1
   
    %%to be removed
    if ~strcmpi(S.type(end-2:end),'dir') 
      S1=spec2spec(S,strcat(S.type(1:end-3),'dir'));
    end
    if ~isfield(S,'phi') || isempty(S.phi), S1.phi=0;
      
    end
    th=S1.theta(:)-S1.phi;
    Sw=simpson(th,S1.S).'; % integral S(w,th) dth
    m=simpson(S1.w,Sw);    % the variance
    mtext={'m0'};
    
    if nr>0
      w=S1.w(:);
      nw=length(w);
      if strcmpi(vari(1),'x')
        Sc=simpson(th,S1.S.*(cos(th)*ones(1,nw))).';
        % integral S*cos(th) dth
      end
      if strcmpi(vari(1),'y')
        Ss=simpson(th,S1.S.*(sin(th)*ones(1,nw))).';
        % integral S*sin(th) dth
        if strcmpi(vari(1),'x')
        Sc=simpson(th,S1.S.*(cos(th)*ones(1,nw))).';
        end
      end
      if ~isfield(S1,'g')
        S1.g=gravity;
      end
      kx=w.^2/S1.g(1); % maybe different normalization in x and y => diff. g
      ky=w.^2/S1.g(end);
      
      if Nv>=1
        switch vari
          case 'x'
            vec = kx.*Sc;
            mtext(end+1)={'mx'};
          case 'y'
            vec = ky.*Ss;
            mtext(end+1)={'my'};
          case 't'
            vec = w.*Sw;
           mtext(end+1)={'mt'}; 
        end 
      else
        vec = [kx.*Sc ky.*Ss w*Sw];
        mtext(end+(1:3))={'mx', 'my', 'mt'};
      end
      if nr>1
      if strcmpi(vari(1),'x')
        Sc=simpson(th,S1.S.*(cos(th)*ones(1,nw))).';
        % integral S*cos(th) dth
        Sc2=simpson(th,S1.S.*(cos(th).^2*ones(1,nw))).';
        % integral S*cos(th)^2 dth
      end
      if strcmpi(vari(1),'y')||strcmpi(vari(2),'y')
        Ss=simpson(th,S1.S.*(sin(th)*ones(1,nw))).';
        % integral S*sin(th) dth
        Ss2=simpson(th,S1.S.*(sin(th).^2*ones(1,nw))).';
        % integral S*sin(th)^2 dth
        if strcmpi(vari(1),'x')
          Scs=simpson(th,S1.S.*((cos(th).*sin(th))*ones(1,nw))).';
          % integral S*cos(th)*sin(th) dth
        end
      end
      if ~isfield(S1,'g')
        S1.g=gravity;
      end
        
      if Nv==2
        switch vari
          case 'xy'
            vec=[kx.*Sc ky.*Ss kx.^2.*Sc2 ky.^2.*Ss2 kx.*ky.*Scs];
            mtext(end+(1:5))={'mx','my','mxx', 'myy', 'mxy'};
          case 'xt'
            vec=[kx.*Sc w.*Sw kx.^2.*Sc2 w.^2.*Sw kx.*w.*Sc];
            mtext(end+(1:5))={'mx','mt','mxx', 'mtt', 'mxt'};
          case 'yt'
            vec=[ky.*Ss w.*Sw ky.^2.*Ss2 w.^2.*Sw ky.*w.*Ss];
            mtext(end+(1:5))={'my','mt','myy', 'mtt', 'myt'};
        end
      else
        vec=[kx.*Sc ky.*Ss w.*Sw kx.^2.*Sc2 ky.^2.*Ss2  w.^2.*Sw kx.*ky.*Scs kx.*w.*Sc ky.*w.*Ss];
        mtext(end+(1:9))={'mx','my','mt','mxx', 'myy', 'mtt', 'mxy', 'mxt', 'myt'};
      end
      if nr>3
        if strcmpi(vari(1),'x')
          Sc3=simpson(th,S1.S.*(cos(th).^3*ones(1,nw))).';
          % integral S*cos(th)^3 dth
          Sc4=simpson(th,S1.S.*(cos(th).^4*ones(1,nw))).';
          % integral S*cos(th)^4 dth
        end
        if strcmpi(vari(1),'y')||strcmpi(vari(2),'y')
          Ss3=simpson(th,S1.S.*(sin(th).^3*ones(1,nw))).';
          % integral S*sin(th)^3 dth
          Ss4=simpson(th,S1.S.*(sin(th).^4*ones(1,nw))).';
          % integral S*sin(th)^4 dth
          if strcmpi(vari(1),'x')  %both x and y
            Sc2s=simpson(th,S1.S.*((cos(th).^2.*sin(th))*ones(1,nw))).';
            % integral S*cos(th)^2*sin(th) dth
            Sc3s=simpson(th,S1.S.*((cos(th).^3.*sin(th))*ones(1,nw))).';
            % integral S*cos(th)^3*sin(th) dth
            Scs2=simpson(th,S1.S.*((cos(th).*sin(th).^2)*ones(1,nw))).';
            % integral S*cos(th)*sin(th)^2 dth
            Scs3=simpson(th,S1.S.*((cos(th).*sin(th).^3)*ones(1,nw))).';
            % integral S*cos(th)*sin(th)^3 dth
            Sc2s2=simpson(th,S1.S.*((cos(th).^2.*sin(th).^2)*ones(1,nw))).';
            % integral S*cos(th)^2*sin(th)^2 dth
          end
        end
        if Nv==2
          switch vari
            case 'xy'
              vec=[vec kx.^4.*Sc4 ky.^4.*Ss4 kx.^3.*ky.*Sc3s ...
                    kx.^2.*ky.^2.*Sc2s2 kx.*ky.^3.*Scs3];
              mtext(end+(1:5))={'mxxxx','myyyy','mxxxy','mxxyy','mxyyy'};
            case 'xt'
              vec=[vec kx.^4.*Sc4 w.^4.*Sw kx.^3.*w.*Sc3 ...
                    kx.^2.*w.^2.*Sc2 kx.*w.^3.*Sc];
              mtext(end+(1:5))={'mxxxx','mtttt','mxxxt','mxxtt','mxttt'};
            case 'yt'
              vec=[vec ky.^4.*Ss4 w.^4.*Sw ky.^3.*w.*Ss3 ...
                    ky.^2.*w.^2.*Ss2 ky.*w.^3.*Ss];
              mtext(end+(1:5))={'myyyy','mtttt','myyyt','myytt','myttt'};
          end
        else
          vec=[vec kx.^4.*Sc4 ky.^4.*Ss4 w.^4.*Sw kx.^3.*ky.*Sc3s ...
               kx.^2.*ky.^2.*Sc2s2 kx.*ky.^3.*Scs3 kx.^3.*w.*Sc3 ...
               kx.^2.*w.^2.*Sc2 kx.*w.^3.*Sc ky.^3.*w.*Ss3 ...
               ky.^2.*w.^2.*Ss2 ky.*w.^3.*Ss kx.^2.*ky.*w.*Sc2s ...
               kx.*ky.^2.*w.*Scs2 kx.*ky.*w.^2.*Scs];
          mtext(end+(1:15))={'mxxxx','myyyy','mtttt','mxxxy','mxxyy',...
          'mxyyy','mxxxt','mxxtt','mxttt','myyyt','myytt','myttt','mxxyt','mxyyt','mxytt'};  

        end % if Nv==2 ... else ...
      end % if nr>3      
      end % if nr>1
      m=[m simpson(w,vec)];
    end % if nr>0
  %  end; %%if Nv==1... else...    to be removed
end % ... else two-dim spectrum 
