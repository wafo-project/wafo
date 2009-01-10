function [sampl,mu1o, mu1oStd, inds]  = cov2csdat(xo,R,cases,method,inds)
% COV2CSDAT generates conditionally simulated values
%
% CALL:  [sample,mu1o, mu1oStd]  = cov2csdat(x,R,cases,method,inds)
%
%   sample = a random sample of the missing values conditioned on the
%            observed (known) data.
%     mu1o = expected values of the missing values conditioned on the
%            observed data.
%  mu1oStd = Standard deviation of mu1o.
%
%       x  = datavector including missing data. 
%             (missing data must be NaN if inds is not given)
%       R  = Auto Covariance Function structure 
%            (NB! must have the same spacing as x)
%    cases = number of cases, i.e., number of columns of sample (default=1)
%   method = 'approximate': (default) condition only on the closest
%                            points. Pros: quite fast
%            'pseudo': Uses the pseudo inverse to calculate conditional
%                      covariance matrix
%            'exact' : doing the exact simulation.
%                     Cons: Slow for large data sets, may not 
%                           return any result due to the covariance 
%                           matrix being singular or nearly singular.
%     inds = indices to spurious or missing data (see x)
%
%  COV2CSDAT generates the missing values from x conditioned on the observed
%  values assuming x comes from a multivariate Gaussian distribution
%  with zero expectation and Auto Covariance function R.
%
% See also reconstruct, rndnormnd, cov2sdat

%Tested on : matlab 5.3, 5.1
% History:
% revised jr 02.07.2000
%  - mnormrnd -> rndnormnd
% revised pab 17.01.2000
%  - updated documentation
% revised pab 12.10.1999
%  
% last modified by Per A. Brodtkorb 17.09.98,31.08.98,27.08.98,16.08.98

% secret methods:
%         'dec1-3': different decomposing algorithm's 
%                   which is only correct for a variables
%                   having the Markov property 
%                   Cons:  3 is not correct at all, but seems to give
%                          a reasonable result 
%                   Pros: 1 is slow, 2 is quite fast and 3 is very fast
%                   Note: (mu1oStd is not given for method ='dec3')

x=xo(:);
N=length(x);
if isstruct(R)
  ACF=R.R;
else
  error('Covariance must be a struct')
end
[n m]=size(ACF);

if (n==m)&&(m~=1),%No Auto Covariance Function is given
  error('not an ACF')
end
ACF=ACF(:);n=max(n,m);
[y I]=max(ACF);
switch I,
  case 1,% ACF starts with zero lag    
  otherwise,error('This is not a valid ACF!!')
end
  
if (nargin==5)&&~isempty(inds),
  x(inds)=NaN;
end
inds=isnan(x);%indices to the unknown observations
Ns=sum(inds);% # missing values
if Ns==0,
 disp('Warning: No missing data, unable to continue.')
 sampl=[];
 return
end
if nargin<3||isempty(cases),
  cases=1;
end
if nargin<4||isempty(method),
  method='approximate'; % default method
end

if Ns==N,% simulated surface from the apriori distribution
  sampl=cov2sdat(R,[N cases]);
  disp('All data missing,  returning sample from the unconditional distribution.')
  return
end

%initializing variables
mu1o=zeros(Ns,1);
mu1oStd=mu1o;
sample=zeros(Ns,cases);
if method(1)=='d',
  xs=cov2sdat(R,[N cases]);% simulated surface from the apriori distribution
  mu1os=sample;
end


switch method(1:4),
  case 'dec1' , % only correct for variables having the Markov property
    % but still seems to give a reasonable answer. Slow procedure.
    Sigma=sptoeplitz([ACF;zeros(N-n,1)]);
 
    %Soo=Sigma(~inds,~inds); % covariance between known observations
    %S11=Sigma(inds,inds); % covariance between unknown observations
    %S1o=Sigma(inds,~inds);% covariance between known and unknown observations
    %tmp=S1o*pinv(full(Soo)); 
    %tmp=S1o/Soo; % this is time consuming if Soo large
    tmp=2*Sigma(inds,~inds)/(Sigma(~inds,~inds) +Sigma(~inds,~inds)' );
    
    if nargout==3, %|nargout==0,
      %standard deviation of the expected surface
      %mu1oStd=sqrt(diag(S11-tmp*S1o'));
      mu1oStd=sqrt(diag(Sigma(inds,inds)-tmp*Sigma(~inds,inds)));
    end
    
    %expected surface conditioned on the known observations from x
    mu1o=tmp*(x(~inds));
    %expected surface conditioned on the known observations from xs
    mu1os=tmp*(xs(~inds,:));
    % sampled surface conditioned on the known observations
    sample=mu1o(:,ones(1,cases))+xs(inds,:)-mu1os; 
    
  case 'dec2',% only correct for variables having the Markov property
    % but still seems to give a reasonable answer
    % approximating the expected surfaces conditioned on 
    % the known observations from x and xs by only using the closest points
    Sigma=sptoeplitz([ACF;zeros(n,1)]);
    n2=floor(n/2);
    inds2=find(inds);
    idx=(1:2*n)'+max(0,inds2(1)-n2);% indices to the points used
    tmpinds=inds; % temporary storage of indices to missing points
    tinds=tmpinds(idx);% indices to the points used
    ns=sum(tinds); % number of missing data in the interval
    nprev=0; % number of previously simulated points
    xsinds=xs(inds,:);
    while ns>0,
      tmp=2*Sigma(tinds,~tinds)/(Sigma(~tinds,~tinds)+Sigma(~tinds,~tinds)');
      if nargout==3||nargout==0,
	%standard deviation of the expected surface
	%mu1oStd=sqrt(diag(S11-tmp*S1o'));
	mu1oStd((nprev+1):(nprev+ns))=...
	    max(mu1oStd((nprev+1):(nprev+ns)), ...
	    sqrt(diag(Sigma(tinds,tinds)-tmp*Sigma(~tinds,tinds))));
      end
      
      %expected surface conditioned on the closest known observations
      % from x and xs2
      mu1o((nprev+1):(nprev+ns))=tmp*(x(idx(~tinds)));
      mu1os((nprev+1):(nprev+ns),:)=tmp*(xs(idx(~tinds),:));      
       
      if idx(end)==N,% 
	ns=0; % no more points to simulate
      else
	% updating by  putting expected surface into x 	
	x(idx(tinds))=mu1o((nprev+1):(nprev+ns));
	xs(idx(tinds))=mu1os((nprev+1):(nprev+ns));

	nw=sum(tinds((end-n2):end));% # data which we want to simulate once 
	tmpinds(idx(1:(end-n2-1)))=0; % removing indices to data ..
	                              % which has been simulated
	nprev=nprev+ns-nw;% update # points simulated so far
				      
	if (nw==0)&&(nprev<Ns), 
	  idx=(1:2*n)'+(inds2(nprev+1)-n2); % move to the next missing data
	else
	  idx=idx+n;
	end
	tmp=N-idx(end);
	if tmp<0, % checking if tmp exceeds the limits
	  idx=idx+tmp;
	end
	% find new interval with missing data
	tinds=tmpinds(idx);
	ns= sum(tinds);% # missing data
      end  
    end
    % sampled surface conditioned on the known observations
    sample=mu1o(:,ones(1,cases))+xsinds-mu1os; 
    
  case 'dec3', % this is not correct for even for variables having the 
    % Markov property but still seems to give a reasonable answer
    % a quasi approach approximating the expected surfaces conditioned on 
    % the known observations from x and xs with a spline
    inds2=find(inds);indg=find(~inds);
    mu1o=interp1(indg, x(~inds),inds2,'spline');
    mu1os=interp1(indg, xs(~inds,:),inds2,'spline');
    % sampled surface conditioned on the known observations
    sample=mu1o(:,ones(1,cases))+xs(inds,:)-mu1os; 

  case {'exac','pseu'}, % exact but slow. It also may not return any result
    Sigma=sptoeplitz([ACF;zeros(N-n,1)]);
    %Soo=Sigma(~inds,~inds); % covariance between known observations
    %S11=Sigma(inds,inds); % covariance between unknown observations
    %S1o=Sigma(inds,~inds);% covariance between known and unknown observations
    %tmp=S1o/Soo; % this is time consuming if Soo large
    if method(1)=='e',%exact
      tmp=2*Sigma(inds,~inds)/(Sigma(~inds,~inds)+Sigma(~inds,~inds)');
    else % approximate the inverse with pseudo inverse
      tmp=Sigma(inds,~inds)*pinv(Sigma(~inds,~inds));
    end
    %expected surface conditioned on the known observations from x
    mu1o=tmp*(x(~inds));
    % Covariance conditioned on the known observations
    Sigma1o=Sigma(inds,inds)-tmp*Sigma(~inds,inds);
    %sample conditioned on the known observations from x
    sample=rndnormnd(mu1o,Sigma1o,cases );
    
    if nargout==3, %|nargout==0,
      %standard deviation of the expected surface
      mu1oStd=sqrt(diag(Sigma1o));
    end
     case 'circ',% approximating by embedding an circulant matrix
       % using that the inverse of the circulant covariance matrix has 
       % approximately the same bandstructure as the inverse of the
       % covariance matrix
       xx2=[x ;zeros(2*n,1);x(end-1:2)];   
       inds=inds(:);
       inds2=[inds ;ones(2*n,1); inds(end-2:2)]; 
       nfft=2^nextpow2(2*N+2*n-2);
       acfC=[ACF(1:n);zeros(nfft-2*n+2,1);ACF((n-1):-1:2)];% circulant vector
       % eigenvalues to the circulant covariance matrix
       lambda=real(fft(acfC,nfft));
       % eigenvalues to the inverse circulant covariance matrix
       lambdainv=1./lambda;
       k1=find(isinf(lambdainv));
       if any(k1),
         disp('Warning procedure is not accurate')
         lambdainv(k)=0;
       end
       acfCinv=real(fft(lambdainv,nfft)/nfft)'; % circulant vector
       k=find(abs(acfCinv/acfCinv(1))>0.2);
       sigmainv = spdiags( acfCinv(ones(nfft,1),k), k-1, nfft, nfft);
      %not finished
  case 'appr',% approximating by only  condition on 
    % the closest points
    % checking approximately how many lags we need in order to 
    % ensure conditional independence
    % using that the inverse of the circulant covariance matrix has 
    % approximately the same bandstructure as the inverse of the
    % covariance matrix
    if 0,
      nfft=2^nextpow2(2*N-2);
      acfC=[ACF(1:n);zeros(nfft-2*n+2,1);ACF((n-1):-1:2)];% circulant vector
      % eigenvalues to the circulant covariance matrix
      lambda=real(fft(acfC,nfft));
      % eigenvalues to the inverse circulant covariance matrix
      lambda=1./lambda;
      lambda(isinf(lambda))=0;
      acfC=real(fft(lambda,nfft)/nfft); % circulant vector
      acfC=acfC/acfC(1); % relative significance of each lag
      % minimum lag distance for which samples are conditionally 
      % independent, Nsig is set to 4 times this value
      Nsig=find(abs(acfC(1:nfft/2))>0.02);       
      Nsig=min(max(2*n,4*Nsig(end)+2),N);% size of Sigma
    else
      Nsig=2*n;
    end
    Sigma=sptoeplitz([ACF;zeros(Nsig-n,1)]);
    n2=floor(Nsig/4);
    inds2=find(inds);
    idx=(1:Nsig)'+max(0,inds2(1)-n2);% indices to the points used
    tmpinds=inds; % temporary storage of indices to missing points
    tinds=tmpinds(idx);% indices to the points used
    ns=sum(tinds); % number of missing data in the interval
    nprev=0; % number of previously simulated points
    x2=x;
    
    while ns>0,
      %make sure MATLAB uses a symmetric matrix solver
      tmp=2*Sigma(tinds,~tinds)/(Sigma(~tinds,~tinds)+Sigma(~tinds,~tinds)');
      Sigma1o=Sigma(tinds,tinds)-tmp*Sigma(~tinds,tinds);
      if nargout==3||nargout==0,
	%standard deviation of the expected surface
	%mu1oStd=sqrt(diag(S11-tmp*S1o'));
	mu1oStd((nprev+1):(nprev+ns))=...
	    max( mu1oStd((nprev+1):(nprev+ns)) , sqrt(diag(Sigma1o)));
      end

      %expected surface conditioned on the closest known observations from x
      mu1o((nprev+1):(nprev+ns))=tmp*(x2(idx(~tinds)));
      %sample conditioned on the known observations from x
      sample((nprev+1):(nprev+ns),:) =...
	  rndnormnd(tmp*(x(idx(~tinds))),Sigma1o,cases );     
      if idx(end)==N,% 
	ns=0; % no more points to simulate
      else
	% updating
	x2(idx(tinds))= mu1o((nprev+1):(nprev+ns)); %expected surface
	x(idx(tinds))=sample((nprev+1):(nprev+ns),1);%sampled surface 
	nw=sum(tinds((end-n2):end));% # data we want to simulate once more
	tmpinds(idx(1:(end-n2-1)))=0; % removing indices to data ..
	                              % which has been simulated
	nprev=nprev+ns-nw;% update # points simulated so far
	
	if (nw==0)&&(nprev<Ns), 
	  idx=(1:Nsig)'+(inds2(nprev+1)-n2); % move to the next missing data
	else
	  idx=idx+n;
	end
	tmp=N-idx(end);
	if tmp<0, % checking if tmp exceeds the limits
	  idx=idx+tmp;
	end
	% find new interval with missing data
	tinds=tmpinds(idx);
	ns= sum(tinds);% # missing data in the interval
      end
    end
end

%size(mu1oStd)
if nargout==0
  plot(find(~inds),x(~inds),'.')
  hold on,
  ind=find(inds);
  plot(ind,mu1o   ,'*')
  plot(ind,sample,'r+')
  %mu1oStd
  plot(ind,[mu1o-2*mu1oStd mu1o+2*mu1oStd ] ,'d')
  %plot(xs),plot(ind,mu1os,'r*')
  hold off
  legend('observed values','mu1o','sampled values','2 stdev')
  %axis([770 850 -1 1])
  %axis([1300 1325 -1 1])
else
  sampl=sample;
end

function y=sptoeplitz(x)
  k=find(x);
  x=x(:)';
  n=length(x);
  if length(k)>0.3*n,
    y=toeplitz(x);
  else
    y = spdiags( x(ones(n,1),k), k-1, n, n);
    k(k(1)==1)=[];
    y=y+spdiags( x(ones(n,1),k), -k+1, n, n);
 end
