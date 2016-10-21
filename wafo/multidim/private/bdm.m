function DS = bdm(Sxyn,Gwt,theta,fi,k,opt)
%BDM  Bayesian Directional Spectrum Estimation Method 
%
%  CALL: DS = bdm(Sxyn,Gwt,thetai,fi,k)
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
%By pab 2004
% based on function BDM in DIWASP toolbox


% Copyright (C) 2000  Per A. Brodtkorb
% 
%  This file, BDM.M, is part of WAFO.
% 
%     BDM is free software; you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation; either version 3 of the License, or
%     (at your option) any later version.
% 
%     BDM is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
    


maxIter = 10;
nmodmin  = 1;
nmodmax  = 3;
Li = 1;
errorTol = 0.001;
if nargin>5
  if ~isempty(opt.errortol), errorTol = opt.errortol; end
  if ~isempty(opt.maxiter),  maxIter  = opt.maxiter; end
  if ~isempty(opt.relax),    Li       = opt.relax; end
  if ~isempty(opt.maxcoef),  maxCoef  = opt.maxcoef; end
  if ~isempty(opt.message),  display  = opt.message; end
  if ~isempty(opt.coefabstol),    coefAbsTol = opt.coefabstol; end
  if ~isempty(opt.minmodelorder), nmodmin = opt.minmodelorder; end
  if ~isempty(opt.maxmodelorder), nmodmax = opt.maxmodelorder; end
end
[m, nt, nf] = size(Gwt);

DS = repmat(1/(2*pi),nt,nf);
H   = zeros(m*m,nt,nf);
phi = zeros(m*m,nf);

% Eliminate meaningless equations such as those determined 
% from the zero co-spectrum and zero quadrature-spectrum.
M = 0;  % M is the number of independent equations

dtheta=theta(2)-theta(1);
tol = sqrt(eps); % treshold defining zero for transfer functions
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

%pause

H   = H(1:M,:,:);
phi = phi(1:M,:);


  
dd            = diag(ones(nt,1))+diag(-2*ones(nt-1,1),-1)+diag(ones(nt-2,1),-2);
dd(1,nt-1:nt) = [1 -2];
dd(2,nt)      = 1;
   
  
h = waitbar(0,'Please wait...BDM calculation');
for ff=k
  waitbar(ff/k(end),h)
  
  A = H(1:M,1:nt,ff)*dtheta;
  B = phi(:,ff);
  n = nmodmin-1;
  stop = 0;
  
  while(~stop)
    n = n+1;
    
    u = 0.5^n; %Hyperparameter minimizing the ABIC (found by trial and error)
    x = repmat(log(1/(2*pi)),nt,1);
    stddiff = 1;
    rlx = Li;
    count = 0;
    while(stddiff>errorTol),
      count            = count+1;
      F                = exp(x);
      E                = diag(F);
      A2               = A*E;
      B2               = B-A*F+A2*x;
      Z(1:M,1:nt)      = A2;
      Z(M+1:M+nt,1:nt) = u*dd;
      Z(1:M,nt+1)      = B2;
      Z(M+1:M+nt,nt+1) = zeros(nt,1);
      
      [Q,UZ] = qr(Z);
      
      TA = UZ(1:nt,1:nt);
      Tb = UZ(1:nt,nt+1);
      
      x1 = TA\Tb;
      stddiff = std(x-x1);
      x = (1-rlx)*x+rlx*x1;
      
      if(count>maxIter || sum(isfinite(x))~=nt)
        if(rlx>Li*2^-4)
          rlx=rlx*0.5;
	  
          if(sum(isfinite(x))~=nt)
            x = repmat(log(1/(2*pi)),nt,1);
          end
          count=0;
        else
          stop = (n>1);
          break;
        end
      end
    end % while(stddiff>0.001)
    
    sig2    = ((norm(A2*x-B2)).^2+(u*norm(dd*x)).^2)/M;
    ABIC(n) = M*(log(2*pi*sig2)+1)-nt*log(u*u)+sum(log(diag(TA).^2));
    
    if(n>1)
      if(ABIC(n)>ABIC(n-1))
        stop = 1;
        n    = n-1;
      end
    end
         
    if(~stop)
      xold=x;
    else
      x=xold;
    end
    
    stop = (n>=nmodmax);
  end % while ~stop
  
  DS(:,ff)=exp(x);
end % for ff = k

close(h)
warning('on');
DS = normspfn(DS,theta);

return


