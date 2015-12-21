function DS = imlm(Sxy,Gwt,thetai,fi,k,opt)
%IMLM  Iterated maximum likelihood method for estimating the directional distribution
%
% CALL  DS = imlm(Sxy,Gwt,thetai,fi,k);
%
%  DS     = Directional distribution (spreading function) size nt x nf
%  Sxy    = matrix of cross spectral densities size m x m x nf
%  Gwt    = matrix of transfer function (abs(Gwt)==1) size m x nt x nf
%  thetai = angle vector length nt
%  fi     = frequency vector length nf
%  k      = index vector to frequencies where Sf>0 length <= nf
%
%  (m  = number of measurement devices)
%  nf  = number frequencies (f or w)
%  nt  = number of angles   (theta)
% 

  
% revised pab jan2005
% Added opt to input  
[m,nt,nf] = size(Gwt);

DS = mlm(Sxy,Gwt,thetai,fi,k);

DS0   = DS;
DSold = DS;
% Parameters controlling the the convergence 
% Li = relaxation parameter (0< Li<=1.5) 1..1.2 is proposed by Krogstad
% Bi = exponent ( 0< Bi) Bi = 1..3 seems appropriate
Li = 1.4; 
Bi = 2;          
errorTol = 1e-3; % maximum tolerance
maxIter   = 30;     % maximum number of iterations
display =0;
if nargin>5
  if ~isempty(opt.errortol), errorTol = opt.errortol; end
  if ~isempty(opt.maxiter),  maxIter  = opt.maxiter; end
  if ~isempty(opt.relax),    Li       = opt.relax; end
  if ~isempty(opt.message),  display  = opt.message; end
end
tolold = inf;
    
h = waitbar(0,'Please wait...IMLM calculation');
for iz = 1:maxIter
  waitbar(iz/maxIter,h)
  % Calculation of cross spectra based on DS
  for ix=1:m
    Sxy(ix,ix,:) = simpson(thetai,squeeze(Gwt(ix,:,:).*conj(Gwt(ix,:,:))).*DS);
    for iy=(ix+1):m,
      Sxy(ix,iy,:) = simpson(thetai,squeeze(Gwt(ix,:,:).*conj(Gwt(iy,:,:))).*DS); 
      Sxy(iy,ix,:) = conj(Sxy(ix,iy,:));
    end
  end
  tmp = (DS0-mlm(Sxy,Gwt,thetai,fi,k));
  DS  = DS+Li*sign(tmp).*abs(tmp.^Bi);
  %Normalize so that int D(theta,f) dtheta = 1 for each f 
  %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  DS   = normspfn(DS,thetai);
  tol  = max(abs(DS(:)-DSold(:)));
  disp(['Iteration nr ' num2str(iz),' of ' num2str(maxIter),' Error = ' num2str(tol) ])
  if iz>5,
    if (min(tol,abs(tol-tolold)*3) < errorTol) , 
      disp('Close enough to convergence'),break,
    end
  end
  tolold = tol;
  DSold  = DS;
    
%   if 0, % used for debugging
%     Di = createspec('dir','f');
%     Di.f = fi;
%     Di.S = DS;
%     Di.theta = thetai;
%     plotspec(Di,2)
%     pause
%   end
end     
close(h)
return; % imlm
  