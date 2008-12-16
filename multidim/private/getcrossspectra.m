function Sxy = getCrossSpectra(thetai,Gwt,DS)
%GETCROSSSPECTRA Compute the cross spectra by integration
%
%  CALL: Sxy = getCrossSpectra(thetai,Gwt,DS);
%
%  thetai = angle vector
%  Gwt    = matrix of transfer functions
%  DS     = directional spectrum
    
  [m,nt,nf] = size(Gwt);
  Sxy = zeros(m,m,nf);
  for ix=1:m
    Sxy(ix,ix,:) = simpson(thetai,squeeze(Gwt(ix,:,:).*conj(Gwt(ix,:,:))).*DS);
    for iy=(ix+1):m,
      Sxy(ix,iy,:) = simpson(thetai,squeeze(Gwt(ix,:,:).*conj(Gwt(iy,:,:))).*DS); 
      Sxy(iy,ix,:) = conj(Sxy(ix,iy,:));
    end
  end
  return