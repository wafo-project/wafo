function [f_mM, paramu,logstr] = wminmax(spec,nit,paramu,t)
%WMINMAX Calculates joint density of minimum and following maximum
%        in a zero-mean stationary Gaussian with normalized spectrum spec
%
% CALL:  f_mM   = wminmax(spec,nit,paramu,t);     
%
%        f_mM  = joint density of minimum and following maximum in  X(t),
%                for def=0 otherwise f_hh is the joint density of (H_1,H_2)
%                that is heights of crest and trough in up-crossing waves.
%                
%        ALL THIS INPUTS NEEDS TO BE SET (NO DEFAULT VALUES ARE ALLOWED)
%    spec   = normalized spectrum L0=L2=1
%    nit    = order of numerical integration: 0,1,2,3,4,5.
%    paramu = parameter vector defining discretization of min/max values.
%    t      = grid of time points between maximum and minimum (to
%             integrate out). interval between maximum and the following
%             minimum,
  
%History
% revised pab Dec2003
% -  replaced code with call to spec2cov2

if paramu(3)<1
  error('Require n>0.')
end
paramv = paramu;
par    = [paramu, paramv];
%  The variable ISQ marks which type of conditioning will be used ISQ=0
%  means random time where the probability is minimum, ISQ=1 is the time
%  where the variance of the residual process is minimal(ISQ=1 is faster).
%
%  NIT, IAC are  described in CROSSPACK paper, EPS0 is the accuracy constant
%  used in choosing the number of nodes in numerical integrations
%  (XX1, H1 vectors). The nodes and weights and other parameters are
%  read in the subroutine INITINTEG from files Z.DAT, H.DAT and ACCUR.DAT.
%
%
%    NIT=0, IAC=1 then one uses RIND0 - subroutine, all other cases
%    goes through RIND1, ...,RIND5. NIT=0, here means explicite formula
%    approximation for XIND=E[Y^+1{ HH<BU(I)<0 for all I, I=1,...,N}], where
%    BU(I) is deterministic function.
%
%    NIT=1, leads tp call RIND1, IAC=0 is also explicit form approximation,
%    while IAC=1 leads to maximum one dimensional integral.
%    .......
%    NIT=5, leads tp call RIND5, IAC is maximally 4-dimensional integral,
%    while IAC=1 leads to maximum 5 dimensional integral.


IAC = 1;
ISQ = 0;
EPS = 5e-5;
EPSS = 1e-6;
EPS0 = 1e-5;
tol = [IAC,ISQ,EPS,EPSS,EPS0];
%tol=[1, 0, 5e-5, 1e-6, 1e-5];

%g=[(-5:0.02:5)', (-5:0.02:5)']; 
 
g = spec.tr;

if length(t)>101
  error('nr. of time points limited to 101.')
end


if abs(t(1))>0.00001 
  error('t(1) < or > 0.')
end
if length(t) < 2
  error('nr. of wavelength <2.')
end



accur = [nit tol];

dt    = min(diff(t)); % dt = t(2)-t(1);
if dt<=0
  error('dt must be positive')
end

Nt = length(t)-1;
nr = 4;
R = spec2cov2(spec,nr,Nt,dt);
ti = linspace(0,Nt*dt,Nt+1).';
filename = writecov(R,ti,nr);

filename2 = {'t.in','transf.in','accur.in','Mm.in'};
cleanup(filename2{:})

disp('   Writing data.')

  fid = fopen('t.in','wt');
  fprintf(fid,'%12.10f\n',t);
  fclose(fid);

  fid = fopen('accur.in','wt');
  fprintf(fid,'%2.0f %2.0f %2.0f\n',accur(1:3));
  fprintf(fid,'%12.10e %12.10e %12.10e\n',accur(4:6));
  fclose(fid);

  fid = fopen('transf.in','wt');
  fprintf(fid,'%12.10e %12.10e \n',g.');
  fclose(fid);

  fid=fopen('Mm.in','wt');
  fprintf(fid,'%12.10e %12.10e %3.0f\n',par(1:3));
  fprintf(fid,'%12.10e %12.10e %3.0f\n',par(4:6));
  fclose(fid);

disp('   Starting Fortran executable.')
%dos([wafoexepath, 'minmax.exe']);
dos([wafoexepath, 'cov2mmpdfreg.exe']);
disp('   Loading data.')
f_mM = loaddata('Maxmin.out');
f_mM = reshape(f_mM(:,3),paramv(3),paramu(3));
f_mM = rot90(f_mM,-2);

try
  logfile = 'Maxmin.log';
  fid     = fopen(logfile,'r');
  logstr  = fread(fid,inf,'uint8=>char');   % Read file
  fclose(fid);
catch
  warning('WAFO:SPEC2MMPDF','Unable to load the log')
end

cleanup('Maxmin.out','Max.out','min.out',logfile,filename{:},filename2{:})