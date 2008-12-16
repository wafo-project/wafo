function DS = emem(Sxyn,Gwt,theta,fi,k,opt)
%EMEM  Extended Maximum Entropy Method 
%
%  CALL: DS = emem(Sxyn,Gwt,thetai,fi,k)
%

%Reference:
% Hashimoto,N. (1997) 
% "Analysis of the directional wave spectra from field data.", Advances in 
% Coastal and Ocean Engineering, Vol.3., pp.103-143
%
% DIWASP, a directional wave spectra toolbox for MATLAB®: User Manual. 
% Research Report WP-1601-DJ (V1.1), Centre for Water Research, University of 
% Western Australia


%Tested on matlab 6.0  
%History : 
%By pab Oct-Nov 2002
% based on function EMEP in DIWASP toolbox

% Copyright (C) 2000  Per A. Brodtkorb
% 
%  This file, EMEM.M, is part of WAFO.
% 
%     EMEM is free software; you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation; either version 3 of the License, or
%     (at your option) any later version.
% 
%     EMEM is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
    

[m, nt, nf] = size(Gwt);


H   = zeros(m*m,nt,nf);
phi = zeros(m*m,nf);

% Eliminate meaningless equations such as those determined 
% from the zero co-spectrum and zero quadrature-spectrum.
M = 0;  % M is the number of independent equations

dtheta = theta(2)-theta(1);
tol    = sqrt(eps); % treshold defining zero for transfer functions
for ix=1:m
  for iy=ix:m
    Htemp  = Gwt(ix,:,:).*conj(Gwt(iy,:,:));                        
    if (any(any(abs(diff(real(Htemp),1))>tol*dtheta)))
      M        = M+1;
      phi(M,:) = real(Sxyn(ix,iy,:));
      H(M,:,:) = real(Htemp); 
    end
    if any(any(abs(diff(imag(Htemp),1))>tol*dtheta) )
      M            = M+1;
      phi(M,:) = imag(Sxyn(ix,iy,:));
      H(M,:,:) = imag(Htemp);    	  
    end
  end
end

% Note: H and Phi here are normalized in another way than described in Hashimoto,N. (1997).
H   = H(1:M,:,:);
phi = phi(1:M,:);

warningState = warning;
warning('off');

M2 = floor(M/2+1);

% Constants controlling the calculation
coefAbsTol    = 0.01;
coefAbsTol2   = 500;
errorTol      = 0.0005;%sqrt(eps);
maxIter       = 25;
minModelOrder = 1;
maxModelOrder = M2;
Li            = 1 ; %relaxation parameter
AICTol        = 0.1;
maxCoef       = 10000;

display =0;
if nargin>5
  if ~isempty(opt.errortol), errorTol = opt.errortol; end
  if ~isempty(opt.maxiter),  maxIter  = opt.maxiter; end
  if ~isempty(opt.relax),    Li       = opt.relax; end
  if ~isempty(opt.maxcoef),  maxCoef  = opt.maxcoef; end
  if ~isempty(opt.message),  display  = opt.message; end
  if ~isempty(opt.coefabstol),       coefAbsTol = opt.coefabstol; end
  if ~isempty(opt.minmodelorder), minModelOrder = opt.minmodelorder; end
  if ~isempty(opt.maxmodelorder), maxModelOrder = opt.maxmodelorder; end
end



cosn  = cos((1:maxModelOrder)'*theta');
sinn  = sin((1:maxModelOrder)'*theta');
cosnt = permute(repmat(cosn,[1 1 M]),[2 3 1]);
sinnt = permute(repmat(sinn,[1 1 M]),[2 3 1]);

AIC       = zeros(1,maxModelOrder);
XY        = zeros(2*maxModelOrder,M);
coef      = zeros(1,2*maxModelOrder); % [a,b]
deltaCoef = zeros(1,2*maxModelOrder); % [deltaA,deltaB]

Nbest = 0; %Best model order

DS = repmat(1/(2*pi),nt,nf); % initialize DS
     
h = waitbar(0,'Please wait...EMEM calculation');

% compute a fast estimate in order to find a good
% starting guess for the Lagrange multipliers 

%DS0   = log(mlm(Sxyn,Gwt,theta,fi,k));
%Oneth = ones(nt,1);

for ff=k, % loop over frequencies where S(f)>0
  waitbar(ff/k(end),h)
    
  Hj     = H(:,:,ff).';  %H(1:M,1:nt) 
  Phij = repmat(phi(1:M,ff).',nt,1); 
    
  stop = 0;
  
  % Assume model order change smoothly over the frequencies
  localMinModelOrder = max(minModelOrder,ceil(Nbest/2));
  
  N    = localMinModelOrder-1;
     
  while(~stop), % loop until the right model order is found
    N         = N+1; %increment modelOrder
    
    % coef0=starting guess for coefs multipliers 
    coef0 = zeros(1,2*N+1);
    %coef0 = [Oneth, cosn(1:N,:).', sinn(1:N,:).'] \ DS0(:,ff); 
    coef(1:2*N) = coef0(2:2*N+1);
    %coef(1:2*N)      = zeros(1,2*N); % [a,b]
    deltaCoef(1:2*N) = zeros(1,2*N); % [deltaA,deltaB]
    
    %coef(2*N-1:2*N)      = zeros(1,2); % [a,b]
    %deltaCoef(2*N-1:2*N) = zeros(1,2); % [deltaA,deltaB]
    
    count      = 0;
    lambda     = Li; %relaxation parameter	 
    modelFound = 0;
    
    exponent  = (coef(1:N)*cosn(1:N,:)+coef(N+1:2*N)*sinn(1:N,:))';
    a0     = -max(exponent); % trick in order to avoid infinities
    Fn     = exp(a0+exponent);
    Dn     = repmat(Fn/(trapz(Fn)*dtheta),1,M); %Fn(theta|f)/norm(Fn)=Dn(theta|f)
    PhiHD  = (Phij-Hj).*Dn;
    Z      = trapz(PhiHD,1)*dtheta; % 1xM
    error1 = max(abs(Z));
    
    while ( ~stop && ~modelFound),%Use Newton-Raphson iteration to find model.
      
      count = count+1;
      
      %XY is the negative jacobian in the Newton iteration.
      for ix=1:N
        %This better than eq. 97 and 98 in Hashimoto (1997).
        %X(ix,1:M)=
        XY(ix,:)   = (Z.*trapz(Dn.*cosnt(:,:,ix),1) - ...
          trapz(PhiHD.*cosnt(:,:,ix),1))*dtheta;
        %Y(ix,1:M)=
        XY(N+ix,:) = (Z.*trapz(Dn.*sinnt(:,:,ix),1) - ...
		      trapz(PhiHD.*sinnt(:,:,ix),1))*dtheta ;
      end
      
      deltaCoefOld     = deltaCoef(1:2*N);
      deltaCoef(1:2*N) = Z/XY(1:2*N,:); % solve eq. 95 by least square
					
      coef(1:2*N) = coef(1:2*N) + lambda*deltaCoef(1:2*N);
      
      if (maxCoef<inf) %This option is not described in Hashimoto,N. (1997).
        %Make sure the coefficients do not diverge to infinity
        k0 = find(abs(coef(1:2*N))>maxCoef);
        if any(k0)
          deltaCoef(k0)=(sign(deltaCoef(k0)).*maxCoef-...
            (coef(k0)-lambda*deltaCoef(k0)))/lambda;
          coef(k0)     = sign(coef(k0)).*maxCoef;
        end
      end
      
      exponent  = (coef(1:N)*cosn(1:N,:)+coef(N+1:2*N)*sinn(1:N,:))';
      a0     = -max(exponent); % trick in order to avoid infinities
      Fn     = exp(a0+exponent);
      Dn     = repmat(Fn/(trapz(Fn)*dtheta),1,M); %Fn(theta|f)/norm(Fn)=Dn(theta|f)
      PhiHD  = (Phij-Hj).*Dn;
      Z      = trapz(PhiHD,1)*dtheta; % 1xM
      
      error2 = error1;
      error1 = max(abs(Z));
      
      if (~any(abs(deltaCoef(1:2*N))>coefAbsTol)  || (error1 <= errorTol));
        modelFound = 1;
      elseif (count>maxIter || ((error2<error1) && ...
          (any(abs(deltaCoef(1:2*N)-deltaCoefOld)>coefAbsTol2) ))),
	
        % Coefficients are diverging and error is increasing
        % solution: under relax the computation or quit
        if(lambda>Li*2^-4),
          lambda      = lambda*0.5;
          count       = 0;
          deltaCoef(1:2*N) = 0;
          % coef(1:2*N)      = 0;
          coef(1:2*N) = coef0(2:2*N+1);
          Dn     = repmat(1/(2*pi),nt,M); %Dn(theta|f)
          PhiHD  = (Phij-Hj).*Dn;
          Z      = trapz(PhiHD,1)*dtheta; % 1xM
          error2 = inf;
          error1 = max(abs(Z));
        else
          stop = 1;
        end
      end
      
      if 0 %modelFound
        Np = 4;
        subplot(Np,1,1), semilogy(abs(Z)+eps,'*'),hline(errorTol), title('error')
        subplot(Np,1,2),plot(deltaCoef(1:2*N),'g*'), title('deltaCoef')
        subplot(Np,1,3), plot(coef(1:2*N),'r*'), title('coef')
        subplot(Np,1,4), plot(theta,Dn)
        drawnow,%
        disp('Hit any key'),pause
      end
    end %while Newton-Raphson
            
    AIC(N) = M*(log(2*pi*var(Z))+1)+4*N+2; %Aikakes Information Criterion
      
    if(N>localMinModelOrder),
      stop = stop | (AIC(N)+AICTol>AIC(N-1)) ;
    end
    if isnan(AIC(N)), stop = 1;  end %this should never happen
    
    if (stop),
      if (N>localMinModelOrder)
        Nbest = N-1;
        coef(1:2*Nbest) = coefOld;
      else
        Nbest = 1;
        coef = zeros(1,2*Nbest);
      end
    else
      Nbest   = N;
      stop    = (N>=maxModelOrder);
      coefOld = coef(1:2*N);
    end
      
  end % while
  
  if display>0
    disp(sprintf('f = %g \t \t Model order = %d \t \t error = %g',fi(ff),Nbest,max(abs(Z))))
  end
  %plot(deltaCoef,'r*'),hold on,plot(coef(1:2*Nbest),'g*'),hold off,pause
  
  exponent = (coef(1:Nbest)*cosn(1:Nbest,:)+coef(Nbest+1:2*Nbest)*sinn(1:Nbest,:)).';
  a0       = -max(exponent); % trick in order to avoid infinities
  DS(:,ff) = exp(a0+exponent);
        
end % for ff=
close(h)
DS = normspfn(DS,theta);

warning(warningState);

if (all(abs(DS(:)-DS(1))<sqrt(eps)))
  warning('No main direction found. Check the estimated spectrum!')
end

return; % EMEM