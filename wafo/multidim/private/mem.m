function DS = mem(Sxyn,Gwt,thetai,fi,k)
%MEM  maximum entropy method for estimating the directional distribution
%
% CALL:  DS = mem(Sxy,Gwt,thetai,fi,k);
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

%Tested on matlab 6.0  
%History : 
%revised pab Oct-Nov 2002  
%By pab 1999?

[m, nt, nf] = size(Gwt);


H   = zeros(m*m,nt,nf);
phi = zeros(m*m,nf);

% Eliminate meaningless equations such as those determined 
% from the zero co-spectrum and zero quadrature-spectrum.
M = 0;  % M is the number of independent equations

dtheta=thetai(2)-thetai(1);
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
%M
H   = H(1:M,:,:);
phi = phi(1:M,:);

warningState = warning;
warning('off');

% Constants controlling the calculation
coefAbsTol    = 0.1;
coefAbsTol2   = 100;

errorTol      = 0.1;%sqrt(eps);
maxIter       = 25;
Li            = 1 ; %relaxation parameter
maxCoef       = 5000;

display =1;

La(1:M)  = zeros(1,M); % Lagrange multiplier
dLa(1:M) = zeros(1,M); % [deltaA,deltaB]
    

DS = repmat(1/(2*pi),nt,nf); % initialize DS
     
h = waitbar(0,'Please wait...MEM calculation');

% compute a fast estimate in order to find a good
% starting guess for the Lagrange multipliers 
DS0 = log(mlm(Sxyn,Gwt,thetai,fi,k));
Oneth = ones(nt,1);

for ff=k, % loop over frequencies where S(f)>0
  waitbar(ff/k(end),h)
    
  Hj   = H(:,:,ff).';  %H(1:M,1:nt) 
  Phij = repmat(phi(1:M,ff).',nt,1); 
  
  stop = 0;
  count    = 0;
  lambda   = Li; %relaxation parameter  
  La0 = -[Oneth Hj] \ DS0(:,ff); % starting guess for the Lagrange multipliers 
  La(1:M) = La0(2:M+1);
  %La(1:M)  = zeros(1,M); % Lagrange multiplier
  
  dLa(1:M) = zeros(1,M); % [deltaA,deltaB]
    
  %coef(2*N-1:2*N)      = zeros(1,2); % [a,b]
  %deltaCoef(2*N-1:2*N) = zeros(1,2); % [deltaA,deltaB]
  
  %size(La),size(Hj)  
  exponent  = -(Hj*(La.'));
  L0     = -max(exponent); % trick in order to avoid infinities
  Fn     = exp(L0+exponent);
  Dn     = repmat(Fn/(trapz(Fn)*dtheta),1,M); %Fn(theta|f)/norm(Fn)=Dn(theta|f)
  PhiHD  = (Phij-Hj).*Dn;
  B      = trapz(PhiHD,1)*dtheta; % 1xM
    
  error1 = max(abs(B));
    
  while(~stop), %%Use Newton-Raphson iteration to find model.
    count = count+1;
           
    %A is the jacobian in the Newton iteration.
    for ix=1:M
      A(ix,:)  = trapz(PhiHD.*repmat(Hj(:,ix),1,M),1)*dtheta;
    end
      
    %A, disp('Hit a key'),pause
    dLaOld(1:M) = La;
    %dLa(1:M)    = B/(A.'); % solve eq. 27 
    dLa(1:M)   = B*pinv(A.'); % solve eq. 27 
    La(1:M) = La + lambda*dLa;
      
    if (maxCoef<inf) 
      %This option is not described in Hashimoto,N. (1997).
      %Make sure the coefficients do not diverge to infinity
      k0 = find(abs(La)>maxCoef);
      if any(k0)
        dLa(k0)=(sign(dLa(k0)).*maxCoef-...
          (La(k0)-lambda*dLa(k0)))/lambda;
        La(k0)     = sign(La(k0)).*maxCoef;
      end
    end
    exponent  = -(Hj*(La.'));
    L0     = -max(exponent); % trick in order to avoid infinities
    Fn     = exp(L0+exponent);
    Dn     = repmat(Fn/(trapz(Fn)*dtheta),1,M); %Fn(theta|f)/norm(Fn)=Dn(theta|f)
    PhiHD  = (Phij-Hj).*Dn;
    B      = trapz(PhiHD,1)*dtheta; % 1xM
    
    
    error2 = error1;
    error1 = max(abs(B));
      
    if (~any(abs(dLa(1:M))>coefAbsTol)  || (error1 <= errorTol));
      stop = 1;
    elseif (count>maxIter || ((error2<error1) && ...
        (any(abs(dLa-dLaOld)>coefAbsTol2) ))),
      dLa(1:M) = 0;
      La(1:M)  = La0(2:M+1);
      
      % Coefficients are diverging and error is increasing
      % solution: under relax the computation or quit
      if(lambda>Li*2^-4),
        lambda   = lambda*0.5;
        count    = 0;
	  
	  
        Dn     = repmat(1/(2*pi),nt,M); %Dn(theta|f)
        PhiHD  = (Phij-Hj).*Dn;
        B      = trapz(PhiHD,1)*dtheta; % 1xM
        error2 = inf;
        error1 = max(abs(B));
      else
        stop = 1;
      end
    end
      
    if 0 %stop,
      N = 4;
      subplot(N,1,1), semilogy(abs(B)+eps,'*'),hline(errorTol), title('error')
      subplot(N,1,2),plot(dLa,'g*'), hline(coefAbsTol*[-1,1]),title('deltaCoef')
      subplot(N,1,3), plot(La,'r*'), title('coef')
      subplot(N,1,4), plot(thetai,Dn), title('D(theta|f)')
      drawnow,%
      disp('Hit any key'),pause
    end
  end %while Newton-Raphson
            
  
  if display>0
    disp(sprintf('f = %g \t \t Error = %g',fi(ff),error1))
  end
  
  exponent  = -(Hj*(La.'));
  L0        = -max(exponent); % trick in order to avoid infinities
  DS(:,ff)  = exp(L0+exponent);
        
end % for ff=
close(h)
DS = normspfn(DS,thetai);

warning(warningState);

if (all(abs(DS(:)-DS(1))<sqrt(eps)))
  warning('No main direction found. Check the estimated spectrum!')
end

return; % MEM


% % old call kept just in case
% [m, nt, nf] = size(Gwt);
% % This method assumes  a directional spreading of the form
% %   D(theta|w) = C*exp( qj*La.') 
% % where C is a  normalization constant, La is vector of the unknown Lagrange
% % multipliers and qj is a cosine/sinus function.
% 
% % compute a fast estimate in order to find a good
% % starting guess for the Lagrange multipliers 
% DS0 = log(mlm(Sxy,Gwt,thetai,fi,k));
% 
% DS = ones(nt,nf)/(2*pi); % If S(f)==0 then set D(theta,f)=1/(2*pi);
% 
% ii = (1:m)'*ones(1,m);
% jj = ones(m,1)*(1:m);
% 
% [indx, indy] = find(ii < jj);   % indices to upper triangular array
% ind = sub2ind([m m],indx,indy); % convert to linear index
% 
% Sxy = permute(Sxy,[3,1,2]);
% Gwt = permute(Gwt,[2 1 3 ]);
% 
% 
% 
% % Set up the non linear equation to solve or to minimize: 
% % int (Pj(ones(nt,1),:)-qj).*repmat(exp(qj*La.'),1,M2) dtheta^2 =0;
% % simpson(thetai, (Pj(ones(nt,1),:)-qj).*repmat(exp(qj*La.'),1,M2)).^2=0;
% eqstr=inline('((simpson(P1,P2.*repmat(exp(P3*transpose(x)),1,P4))))',4);
% 
% M2 = m*(m-1); M1 = M2/2;
% 
% qj = zeros(nt,M2);
% La = rand(1,M2)*10-5; %zeros(1,M2); % starting guess
% Pj = La;
% if 0
%   tmp = thetai(:) * (1:M1);
%   qj2 = [cos(tmp) sin(tmp)];
%   size(qj2)
% end
% 
% 
% % find matlab version
% vstr = version; vind = find(vstr=='.');
% vstr = str2num(vstr(1:vind(2)-1));
% 
% Oneth = ones(nt,1);
% 
% 
% for ix = k, % looping over non-zero values of S(f) only
%   Pj(1:M1)       = real(Sxy(ix,ind));
%   Pj(M1+1:M2)    = imag(Sxy(ix,ind));
%   tmp            = Gwt(:,indx,ix).*conj(Gwt(:,indy,ix));
%   qj(:,1:M1)     = real(tmp);
%   qj(:,M1+1:M2)  = imag(tmp);
%   %qj2 = qj;
%   %plot(qj)
%   %pause
%   % Find a good starting guess for La based on log(DS_mlm).
%   if 1,
%      La0 = [Oneth qj] \ DS0(:,ix);
%   else
%     La0 = [Oneth qj(:,[1:2 M1+1:M1+2] )] \ DS0(:,ix);
%     La  = [La0(2:3).' zeros(1,M1-2) La0(4:5).' zeros(1,M1-2) ]
%   end
%   %mans = memfun(La,thetai,Pj(ones(nt,1),:)-qj,qj,M2);
%   %if any(isnan(mans)),  La = rand(1,M2)*10-5;  end
%   
%   % find the Lagrange multipliers La 
%  
%   %La = fsolve('memfun',La,optimset('fsolve'),thetai,Pj(ones(nt,1),:)-qj,qj,M2);
%   La = fminsearch('memfun',La,[],thetai,Pj(ones(nt,1),:)-qj,qj(:,1:M2),M2);
%   %La = fminunc('memfun',La,[],thetai,Pj(ones(nt,1),:)-qj,qj,M2)
%  
%   
%   La(isnan(La)) = 0; % make sure that isnans is zero
%   %mans2 = memfun(La,thetai,Pj(ones(nt,1),:)-qj,qj,M2);
%   
%  
%   disp([' Finished ' num2str(ix) ' of ' num2str(nf)])
%   %disp(num2str([mans;mans2]))
%   %disp(num2str(La))
%   DS(:,ix) = exp(qj*La.');
%   %if sum(abs(mans2)) > 1  La = rand(1,M2)*10-5;  end 
% end
% %Normalize so that int D(theta,f) dtheta = 1 for each f 
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% DS = normspfn(DS,thetai);
% 
% return;
