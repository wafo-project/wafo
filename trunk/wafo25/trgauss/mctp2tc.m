function F=mctp2tc(freqPVR,utc,param,freqPVL)
%MCTP2TC  Calculates frequencies for the  upcrossing troughs and crests
%  using Markov chain of turning points.
%
%  CALL: f_tc = mctp2tc(f_Mm,utc,param);
%
%  where
%
%        f_tc  = the matrix with frequences of upcrossing troughs and crests,
%        f_Mm  = the frequency matrix for the Max2min cycles,
%        utc   = the reference level,
%        param = a vector defining the discretization used to compute f_Mm,
%                note that f_mM has to be computed on the same grid as f_mM. 
%
%  optional call: f_tc = mctp2tc(f_Mm,utc,param,f_mM);
%
%        f_mM  = the frequency matrix for the min2Max cycles.
%

error(nargchk(3,4,nargin))
if nargin<4
  freqPVL=freqPVR;
end

u=levels(param);
udisc=fliplr(u);
ntc=sum(udisc>=utc);
%ntc1=sum(u<=utc);

%if udisc(ntc) > utc
%   ntc=ntc+1
%end

%if u(ntc1) >= utc
%   ntc1=ntc1-1
%end


if (ntc < 2) || (ntc > (param(3)-1))
  error('the reference level out of range, stop')
end

% normalization of frequency matrices

nP=length(freqPVR);
for i=1:nP,
  rowsum=sum(freqPVR(i,:));
  if rowsum~=0
    freqPVR(i,:)=freqPVR(i,:)/rowsum;
  end
end
P=fliplr(freqPVR); 

Ph=rot90(fliplr(freqPVL),-1);
for i=1:nP,
    rowsum=sum(Ph(i,:));
    if rowsum~=0
       Ph(i,:)=Ph(i,:)/rowsum;
    end
end
Ph=fliplr(Ph);


n=nP; F=zeros(n,n); 

if ntc > n-1
  error('index for mean-level out of range, stop')
end

 
%F(1:ntc,1:ntc1)=freqPVL(1:ntc,1:ntc1);
F(1:ntc-1,1:(n-ntc))=freqPVL(1:ntc-1,1:(n-ntc));


F=cmat2nt(F);

for i=2:ntc,
	for j=ntc:n-1,

    if i<ntc
  
      Ap=P(i:ntc-1,i+1:ntc); Bp=Ph(i+1:ntc,i:ntc-1);
      dim_p=ntc-i;
      tempp=zeros(dim_p,1);
      I=eye(size(Ap));
      if i==2
        e=Ph(i+1:ntc,1);
      else
        e=sum(Ph(i+1:ntc,1:i-1),2);
      end
      if max(abs(e))>1e-10
        if dim_p==1
          tempp(1,1)=(Ap/(1-Bp*Ap)*e);
        else
          tempp=Ap*((I-Bp*Ap)\e);
	      end
      end
    end

    if j>ntc

      Am=P(ntc:j-1,ntc+1:j); Bm=Ph(ntc+1:j,ntc:j-1);
      dim_m=j-ntc;
      tempm=zeros(dim_m,1);
      Im=eye(size(Am));
      if j==n-1
        em=P(ntc:j-1,n);
      else
        em=sum(P(ntc:j-1,j+1:n),2);
      end
      if max(abs(em))>1e-10
        if dim_m==1
          tempm(1,1)=(Bm/(1-Am*Bm)*em);
        else
          tempm=Bm*((Im-Am*Bm)\em);
        end
      end
    end
     
    if (j>ntc) && (i<ntc)
      F(i,n-j+1)=F(i,n-j+1)+tempp'*freqPVL(i:ntc-1,n-ntc:-1:n-j+1)*tempm;
      F(i,n-j+1)=F(i,n-j+1)+tempp'*freqPVL(i:ntc-1,n-j:-1:1)*ones(n-j,1);
      F(i,n-j+1)=F(i,n-j+1)+ones(1,i-1)*freqPVL(1:i-1,n-ntc:-1:n-j+1)*tempm;
    end
    if (j==ntc) && (i<ntc)
      F(i,n-j+1)=F(i,n-j+1)+tempp'*freqPVL(i:ntc-1,n-j:-1:1)*ones(n-j,1);
      for k=1:ntc
        F(i,n-k+1)=F(i,n-ntc+1);
      end
    end
    if (j>ntc) && (i==ntc)
      F(i,n-j+1)=F(i,n-j+1)+ones(1,i-1)*freqPVL(1:i-1,n-ntc:-1:n-j+1)*tempm;
      for k=ntc:n
        F(k,n-j+1)=F(ntc,n-j+1);
      end
    end
  end
end

F=nt2cmat(F);


%fmax=max(max(F));

%  contour (u,u,flipud(F),...
%fmax*[0.005 0.01 0.02 0.05 0.1 0.2 0.4 0.6 0.8])
%  axis([param(1) param(2) param(1) param(2)])

%  title('Crest-trough density')
%  ylabel('crest'), xlabel('trough')                    
%  axis('square')
%if mlver>1, commers, end








