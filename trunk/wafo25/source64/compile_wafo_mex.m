%% COMPILE_WAFO_MEX Script that compiles WAFO mexfiles.
%
% COMPILE_WAFO_MEX compiles WAFO mexfiles and copies them to the correct
% location. In order to use this you must have a working c and fortran
% compiler and wafo installed.
%
% Example
% cd(fullfile(waforoot,'source'))
% compile_wafo_mex
% 


%localdir = fullfile(waforoot,'source');

%p0 = pwd;
%cd(localdir)

%ext = '.dll'; % for PC compatible with old versions of matlab
ext = ['.' mexext];

outdir = wafomexpath;
if ~exist(outdir,'dir')
   [success,msg,msgid] = mkdir(outdir);
   if ~success
     error(msg)
   end
end

root = waforoot;

answer = input('Compile c-mex files? (y/n) ','s');
if strncmpi(answer,'y',1)
%% C MEX FILE COMPILATION
%The mex files written in c is installed by the following commands in the matlab prompt:
  try
    disp('Choose a c-compiler!')
    mex('-setup')

    cmexfile = {'findcross.c','findrfc.c','disufq.c'};
    coutfile = strrep(cmexfile,'.c',ext);
    
  
    for ix=1:numel(cmexfile)
      mex('-O','-outdir',outdir,'-output',coutfile{ix},cmexfile{ix})
    end
    %clean-mod:
    delete('*.mod')
    disp('Successfully compiled c-mex-files!')
  catch
    disp('Unexpected termination of compilation of c-mex-files!')
  end
end

answer = input('Compile f-mex files? (y/n) ','s');
if strncmpi(answer,'y',1)
%% F MEX FILE COMPILATION
%The mex files written in F is installed by the following commands in the
%matlab prompt:
  try
    disp('Choose a fortran-compiler!')
    mex('-setup')
pwd
  
    fmexfile = {'mvnprodcorrprbmex.F','mvnprdmex.F','mexmvtnormprb.F','mexmvnprb2.F',...
      'mexmvnprb.F','mexGenzMvnPrb.F'};
    foutfile = strrep(fmexfile,'.F',ext);
    
    for ix=1:numel(fmexfile)
      mex('-O','-outdir',outdir,'-output',foutfile{ix},fmexfile{ix})
    end

    fmexfile2 = {'intmodule.f','jacobmod.f', 'rind71.f', 'mexrind71.F';...
      'intmodule.f','jacobmod.f', 'rind2007.f','mexrind2007.F'}; %...
%      'intmodule.f','jacobmod.f', 'alanpab22.f', 'mexrindalan22.f'};
    foutfile2 = strrep(fmexfile2(:,end),'.F',ext);

    for ix=1:size(fmexfile2,1)
      mex('-O','-outdir',outdir,'-output',foutfile2{ix},fmexfile2{ix,:});
    end


%    fmexfile3 = {'mvdist.f90', 'mexmvtprb.f'};
%    foutfile3 = strrep(fmexfile3{end},'.f',ext);
%    mex('-O', '-outdir',outdir,'-output',foutfile3 ,fmexfile3{:});
    
    
    %clean-mod:
    delete('*.mod')
    disp('Successfully compiled f-mex-files!')
  catch
    disp('Unexpected termination of compilation of f-mex-files!')
  end
end
%cd(p0) % move to the starting directory
