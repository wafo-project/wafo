function DS = mlm(Sxy,Gwt,thetai,fi,k,opt)
% MLM  maximum likelihood method for estimating the directional distribution
%
% CALL  DS = mlm(Sxy,Gwt,thetai,fi,k);
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

[m nt nf] = size(Gwt);

 
% size(Sxy),nf,nt
%-----------------------------------------------------
%   inverting matrix of cross-spectra at every frequency
%------------------------------------------------------

I=eye(m)*sqrt(eps);
if 1, % New call, slightly faster 
  % initialize DS
  DS = repmat(1/(2*pi),nt,nf); % If S(f)==0 then set D(theta,f)=1/(2*pi);
  for ix = k(:)', % looping over non-zero values of S(f) only
    %Gmat = Gwt(1:m,1:nt,ix).'; % = transpose(Gwt(1:m,1:nt,ix))
    %H    =  real((conj(Gmat)*pinv(Sxy(:,:,ix)+I)).*Gmat);
     
    H   = real((ctranspose(Gwt(1:m,1:nt,ix))*pinv(Sxy(:,:,ix)+I)).*(Gwt(1:m,1:nt,ix)).');
    tmp = sum(H,2);
    if any(tmp==0)
      if all(tmp==0)
	DS(:,ix)= 1/(2*pi*nt);
      else
	tmp(tmp==0)=eps;
	DS(:,ix) = 1./tmp;
      end
    else
      DS(:,ix) = 1./tmp; %  size(H) = nt x m
    end
  end
else % Old call:
    
  DS  = zeros(nt,nf); %initialize
  for ix=k, % 1:nf
    Sxy(:,:,ix) = pinv(Sxy(:,:,ix)); %Gmn^-1
  end
  for ix=1:m, %    m-1,
    Sm1 = real(squeeze(Sxy(ix,ix,:))).';
    Gm1 = conj(squeeze(Gwt(ix,:,:)));
    DS = DS+Sm1(ones(nt,1),:).*abs(Gm1).^2 ;
    for iy = (ix+1):m, %    m-1,
      Sm1 = squeeze(Sxy(ix,iy,:)).';
      Gm2 = (squeeze(Gwt(iy,:,:)));
      DS  = DS+2*real(Sm1(ones(nt,1),:).*Gm1.*Gm2 );
    end
  end
  DS  = 1./real(DS);
end


% if 0,
%   Di   = createspec('dir','f');
%   Di.f = fi;
%   Di.S = DS;
%   Di.theta = thetai;
%   plotspec(Di,2)
%   pause
% end

%Normalize so that int D(theta,f) dtheta = 1 for each f 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DS = normspfn(DS,thetai);

% if 0,
%   Di = createspec('dir','f');
%   Di.f = fi;
%   Di.S = DS;
%   Di.theta=thetai;
%   plotspec(Di,2)
%   pause
% end
return
