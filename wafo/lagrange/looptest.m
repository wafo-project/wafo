function Nloops = looptest(S,opt,varargin)
%LOOPTEST Simulates 2D Lagrange waves to estimate folding rate
%

Nloops=0;
Nsim=100;
options=simoptset(opt,varargin{:});
for n=1:Nsim,
    waitbar(n/Nsim)
    [~,x]=spec2ldat(S,options);
    if x.folded,
        Nloops=Nloops+1;
    end
end
disp([num2str(Nloops) ' out of ' num2str(Nsim) ' x-fields were folded'])

