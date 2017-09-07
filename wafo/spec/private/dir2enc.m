function Snew=dir2enc(S,d,v)
%DIR2ENC Transform a dir. spectrum to an encountered. (Used in spec2spec)
%
% CALL:  Snew=dir2enc(S,d,v)
%
%      Snew = encountered spectrum, directional or frequency type
%      S    = directional spectrum
%      d    = dimension of result, 1=frequency type, 2=directional type
%                                                 (default 1)
%      v    = speed of ship     (default 0)
%
% Evaluates the spectrum of the seaway encountered by a ship travelling 
% along the x-axis (rotated by in a constant angle (S.phi)) 
% with constant speed (v) from the directional spectrum of the
% seaway.
%

% References: Podgorski et al. (2000)
%             Exact distributions for apparent waves in irregular seas
%             Ocean Engng. 27, p. 979-1016.

% Tested on MATLAB 5.3  
% History: 
% revised:
% pab 21.09.2004: Added todo statements  
%  ir 26.06.04 BUG's removals plus major changes.
%  ir 04.04.01 BUG's removals plus major changes.
%  es 05.06.00 simpson in stead of trapz, clean up 
%  by es 18.08.1999

% TODO % The numerical accuracy of the method is pure. 
% TODO % Should remove singularity by e.g. splitting the integral.

%error(nargchk(1,3,nargin));
narginchk(1,3)
if nargin<0 || isempty(S)
  error('Needs an input spectrum');
end

if ~any(strcmpi(S.type,{'dir','freq'}))
  error('Spectrum must be of type dir or freq ');
end

if (nargin<2 || isempty(d))
  d=1;
end
if ( nargin<3 || isempty(v))
  v=0;
end

Snew=S;
if isfield(Snew,'phi')
  phi=Snew.phi;
else
  phi=0;
  Snew.phi=phi;
end

if ~isfield(Snew,'theta')
  Snew.theta=0;
end

Snew.v=v;
if v<0
  disp('Message: Negative speed, interpreted as opposite direction')
  v=abs(v);
  phi=phi+pi;
end


if (v==0)
  if d==1
    Snew=spec2spec(Snew,'freq');
    Snew.type='enc';
  else
    Snew.type='encdir';
  end
  return
end


th=Snew.theta;
if isfield(Snew,'f')
  w=2*pi*Snew.f(:);
  Sf=Snew.S/2/pi;
else
  w=Snew.w(:);
  Sf=Snew.S;
end

if length(th)>1
  % IR changes 26 june 2004, formula (6) in Podgorski et al. 
  [w1,th1] = meshgrid(w,th);
  Sfneg    = interp2(w,S.theta,Sf,w1,mod(th1+2*pi,2*pi)-pi);
  Sf       = [fliplr(Sfneg(:,2:end)) Sf]/2; % unitary spectum - symmetric in w
  w2       = [-flipud(w(2:end));w];

  % Now we need to make the transformation that 
  %        Sf1(theta,omega)<= Sf(theta+phi,omega).

  [w1,th1] = meshgrid(w2,th);
  
  Sf1      = interp2(w2,S.theta,Sf,w1,mod(th1+phi+pi,2*pi)-pi);

  Snew.S = Sf1;
  g      = gravity;
  % Since the encountered spectrum has usually higher frequencies we
  % need to use longer vector of frequencies omega.

  w        = linspace(0,2*max(w),2*length(w));
  thw      = pi*ones(size(w));
  thw(w>0) = acos(max(-1,-g/4/v./w(w>0))); % This computes integration 
  % limits in  theta
  %plot(w,thw)

  in=zeros(length(th),length(thw));
  for j=1:length(w)                     
    % The limits depends on omega
    % and hence for each omega we
    % check which theta satisfies
    % the integration limits.
    in(abs(th)<thw(j),j)=1;
  end

  in = logical(in);
  s1 = zeros(length(th),length(w2));
  s2 = s1;
  for j=1:length(th)
    ix = find(in(j,:));
    if any(ix)
      h = g/(2*v*cos(th(j)));
      a = sqrt(1+2*w(ix)/h);
      s1(j,ix) = (interp1(w2,Sf1(j,:),h*(-1+a)))./a;
      s2(j,ix) = (interp1(w2,Sf1(j,:),h*(-1-a)))./a;
    end
  end
  s1(isnan(s1)) = 0;
  s2(isnan(s2)) = 0;
  Snew.S        = 2*(s1+s2);
  Snew.w        = w;
  Snew.type     = 'encdir';
  Snew.phi      = 0.;

  if d==1
    Snew.S=simpson(Snew.theta,Snew.S,1); %integrate out angle
  end

else 
  if abs(phi-pi/2)>0.00001
    % method = 'spline';
    method = 'pchip';
    g  = gravity;
    c  = 4*v/g;
    w  = linspace(0,2*max(S.w),2*length(S.w));
    Dpl   = zeros(size(w));
    Dmin  = Dpl;
    %kapa  = phi;      
    h1    = 2./c./cos(phi);
    delta = 1+c*w*cos(phi);
    ind   = find(delta>=0);
    y1    = sqrt(delta(ind));
    lambda1 = h1.*(-1+y1);
    lambda2 = h1.*(-1-y1);
    D_1   = interp1(S.w,S.S,lambda1,method,0)./y1;
    D_2   = interp1(S.w,S.S,lambda2,method,0)./y1; 
    Dpl(ind) = D_1+D_2; 
    
    h1    = 2./c./cos(phi);
    delta = 1-c*w*cos(phi);
    ind   = find(delta>=0);
    y1    = sqrt(delta(ind));
    lambda1 = h1.*(-1+y1);
    lambda2 = h1.*(-1-y1);
    D_1 = interp1(S.w,S.S,lambda1,method,0)./y1;
    D_2 = interp1(S.w,S.S,lambda2,method,0)./y1; 
    %
    % Note that actually the encountered spectrum is always a directional
    % spectrum. 
    %
    Dmin(ind) = D_1+D_2;
    Snew.w    = [fliplr(-w) w];
    Snew.S    = [fliplr(Dmin) Dpl];
    Snew.type = 'encdir';
    if d==1
       Snew.type = 'enc';
       Snew.S    = (Dmin+Dpl)';
       Snew.w    = w';
       if abs(phi)>pi/4
          [mm,im] = max(Snew.S);
          dw0     = abs(S.w(end)-S.w(end-1));
          dw      = abs(w(end)-w(end-1));
          L0      = dw0*sum(S.S);
          %L1=dw*sum(Snew.S)-0.5*mm*dw
          Snew.S(im) = 0.;
          L2         = dw*sum(Snew.S);
          %[L0 L1 2*(L0-L2)/dw]
          Snew.S(im)=max(0,2*(L0-L2)/dw);
       end
    end
  end
  
end

if d==1
  Snew = rmfield(Snew,'theta');
  Snew.type = 'enc';
end
