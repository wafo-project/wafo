function ml2html2(rdir,htmldir,sdef,intyp)
%ml2html2  Produce html files for each m-file in root directory.
%
%  CALL: ml2html2(rdir,htmldir,sdef,typ);
%
% rdir       = Root directory to search for m-files (default pwd)
% htmldir    = Root directory for HTML files. 
%              (default fullfile(pwd,'html','') 
% sdef       = integer defining where to search for m-files:
%              0 Root directory
%              1 Root directory and all subdirectories
%              2 all subdirectories of the root directory
% typ        = [fc hh fl fp fs tt] vector of integer values indicating
%               which part of the m-file should be included in the
%               HTML's. 
%               fc  = 1 if function call (synopsis) is included
%               hh  = 1 if help header is included
%               fl  = 1 if reference to loaded files
%               fp  = 1 if reference to printed files
%               fs  = 1 if reference to saved files
%               tt  = 1 if original m-code is included
%
%ML2HTML2 reads a list of matlab .m files from the directories RDIR
%to produce a hypertext documentation suitable for
%browsing with a WWW browser such as Mosaic or Netscape.  An file
%html_dir/index.html in written. Subdirectories are written for matlab
%directory that is encountered containing
%.html files corresponding to every .m file.
%
% SkipExample
% rdir = fullfile(waforoot, 'fileutil');
% htmldir = [];  
% sdef = 1;
% typ = [1 1 0 0 0]; 
% ml2html2(rdir,htmldir,sdef,typ);
% 
% See also: finddep  

global HTINI
htinit; % Read htini.dat and initialize the toolbox
currdir = pwd;
if nargin<1|isempty(rdir),    rdir    = currdir; end
if nargin<2|isempty(htmldir), htmldir = fullfile(currdir,'html',''); end
if nargin<3|isempty(sdef),    sdef    = 0;end

nt2 = 6; typ = [1 1 0 0 0 0];
if nargin>=4&~isempty(intyp),   
  nt  = length(intyp); 
  ind = find(~isnan(intyp(1:min(nt,nt2))));
  if any(ind),  typ(ind)=intyp(ind);end
end


switch sdef
  case 0, % Find files in root directory only
    udir = rdir;
  case 1, % Root directory and all subdirectories
    udir = findsubdir(rdir); 
    if isempty(udir),return,end
  case 2, % Only subdirectories of the root directory
    udir = findsubdir(rdir);
    if  isempty(udir),return,end
    udir(1:size(rdir,1),:) = []; % remove rdir from list      
end
if isempty(udir),disp('Unable to create HTML files'),return,end

%######################################################################
%#
   disp('Read the m-files, and compute cross-reference information')
%#
%######################################################################

global C W M
if 1,
  [C,W,M] = finddep(udir);
end
C1 = C{1};
C2 = C{2};
[n,m] = size(C1);
M = cellstr(M);
M         = strrep(M,'.M','.m'); % Make sure the ending is lower case
htmlfiles = char(strrep(M,'.m','.htm'));
mfiles    = char(strrep(M,'.m','')); 
M = char(M);

% Removing self-references:
C1 = C1 -speye(n);
C2 = C2 -speye(n);

%#####################################################################
%#
   disp('Setup the html directories')
%#
%######################################################################


% Create the main HTML directory
dirchk = exist(htmldir);
if dirchk ~=7,
  [t, r]= strltok(htmldir,filesep);
  if isempty(r),
    [status msg]=mkdir(htmldir);
  else
    [status msg]=mkdir(r,t);
  end
  if status==0, error(msg),end
  disp(['You did not have the directory ', htmldir ,' so I created it.']);
end;

% copy gif-fil to html directory.
copyfile(fullfile(HTINI.gif_home,HTINI.cp),htmldir); 


% Remove directories with no .m-files
Nw       = length(W);
for ix = Nw:-1:1,
  if isempty(W(ix).m), W(ix)=[];end
end
% Create a html subdirectory name for every unique matlab directory in the
% list W. The name is constructed using the tail of the directory
% prefaced by a unique number.
Nw       = length(W);
[T,R] = strltok(rdir,filesep)
Nz       = max(length(R),1);
htmlsdir = cell(1,Nw);
mdir     = cell(1,Nw);
drmask   = []; % index specifying which directory the HTML files are
for ix = 1:Nw,
  nr       = num2str(ix);
  drmask   = [drmask, ix(ones(1,length(W(ix).m)))];
  mdir(ix) = {W(ix).path(Nz:end)};
  htmlsdir(ix) = {[nr '.' strltok(W(ix).path,filesep)]};
  [status, msg] = mkdir(htmldir,htmlsdir{ix});
  if status==0, error(msg),end
end


%######################################################################
%#
   disp('Write the master index file')
%#
%######################################################################

indexfile = fullfile(htmldir,'index.htm');

disp('Writing master indexfile');

[fid, msg]=fopen(indexfile,'w');
if fid==-1,error(msg),end
fprintf(fid,'%s \n','<TITLE>Matlab Index</TITLE>');
fprintf(fid,'%s \n','<BODY>');
fprintf(fid,'%s \n','<H1>Matlab Index</H1>');
%sigwrt(fid);
%&tagline;

% Print a short introduction

% Print directory listing

fprintf(fid,'%s \n','<HR><H2>Matlab Directory Indices</H2>');
fprintf(fid,'%s \n','<pre> <UL>');


[tmp,ind]=lexsort(char(mdir));
for ix =ind(:).'
  %dr = strltok(W(ix).path,filesep);
  refwrt(fid,mdir{ix},fullfile(htmlsdir{ix},'index.htm') ,1);
end
fprintf(fid,'%s \n','</UL>');

% Include links to every file that was found

fprintf(fid,'%s \n','<HR><H2>Identifiers found in these directories</H2>');

space = ' ';
% We'll do this five across in alphabetical order
[tmp ind] = lexsort(M);
iy=0;
for ix=ind.'
  iy = iy+1; 
  iz = drmask(ix);
  %nr = num2str(iz);
  %dr = strltok(W(iz).path,filesep);
  %filestr = strtok(M(ix,:),'.');
  fileref = fullfile(htmlsdir{iz}, htmlfiles(ix,:)); %[filestr '.htm']);
  %b = space(ones(1,15-length(deblank(mfiles(ix,:)))));
  fprintf(fid,' <A HREF ="%s">%s</A> ',fileref,mfiles(ix,:));
  if mod(iy,5)==0, fprintf(fid,'\n');,end
end;
fprintf(fid,'%s \n','<HR></BODY>');
fclose(fid);

%######################################################################
%#
   disp('Write an index for each html subdirectory')
%#
%######################################################################

for ix=1:Nw
%  nr = num2str(ix);
%  dr = strltok(W(ix).path,filesep);

  [status msg]=mkdir(htmldir,htmlsdir{ix});
  if status==0,error(msg),end
  indexfile = fullfile(htmldir,htmlsdir{ix}, 'index.htm');

  disp(['Writing index file ' indexfile])

  [fid,msg] = fopen(indexfile,'w');
  if fid==-1,error(msg),end
  fprintf(fid,'<TITLE>Index for Directory %s</TITLE>\n',mdir{ix} );
  fprintf(fid,'<BODY>\n');
  fprintf(fid,'<A HREF = \"../index.htm\">[Return to Master Index]</A>\n');
  fprintf(fid,'<H1>Index for %s</H1>\n',mdir{ix});
  sigwrt(fid);

  % Now look for a Readme.m file, seemingly a Matlab standard. If there
  % is one, then the help portion is included in the index file.
  ind = strmatch('readme.m',lower(W(ix).m),'exact');
  for iy=ind.' 
    HH = help(fullfile(W(ix).path,W(ix).m{iy}));
    txtwrt(fid,HH);
  end 

  % Now write the index catalog for the .m files in this directory

  fprintf(fid,'<HR><H2>Matlab files in this Directory</H2>\n<pre>\n');
  iy=0;
  for iz=find(drmask == ix) %1:length(W(ix).m),
    iy = iy+1; 
    %filestr= strtok(W(ix).m{iz},'.');
    %fileref = [filestr '.htm'];
    
    %b = space(ones(1,15-length(deblank(mfiles(iz,:)))));
    fprintf(fid,' <A HREF ="%s">%s</A> ',deblank(htmlfiles(iz,:)),mfiles(iz,:));
    if mod(iy,5)==0, fprintf(fid,'\n');,end
  end;
  fprintf(fid,'</pre>\n');
  fprintf(fid,'<HR></BODY>\n');
  fclose(fid);
end

%######################################################################
%#
      disp('Write an html file for every m-file')
%#
%######################################################################

%# Now write the html file for each matlab file. Need to reread each matlab
%# file to find the help text. Note that we can't do this in a single loop
%# because we want the back reference information, and also some people
%# put the help text before the function declarations.  %

%# Need a list of mfiles with unique identifiers


%%%%%%%%%%%%%  Going into each html file %%%%%%%%%%%%%%%%%
%for i=1:length(drmask),
%  ix = drmask(i);
i = 0;
for ix=1:Nw,
  hdir = strltok(W(ix).path,filesep);
  nr   = num2str(ix);
  for iy =1:length(W(ix).m)
    i=i+1;
    CA = cell(1,6);
    HH = [];FL=[];FP=[];FS=[];tt=[];funcall=[];
    %filename = W(ix).m{iy};
    filename = deblank(M(i,:));
    fullfilename = fullfile(W(ix).path,filename);
    disp(['Check nr:  ' num2str(i) '   ' fullfilename])
    % Creating the string matrices containing
    % names of m-files called by current function or m-files that
    % calls the current function:
    
    cref = cell(1,4); % cellarray of strings with crossreferences
    htmlind = cell(1,4); % cellarray of indices to html directories
    
    % Find functions called by "filename"
    ind = find(C1(i,:)==1); 
    if any(ind), 
      [cref{1} ind0] = lexsort(M(ind,:)); 
      htmlind{1} = drmask(ind(ind0));
    end % ref_to
    
    % Find functions dependent on "filename"
    ind = find(C1(:,i)==1);
    if any(ind), 
      [cref{2},ind0] = lexsort(M(ind,:)); 
      htmlind{2} = drmask(ind(ind0));
    end % 
    
    % Find functions mentioned in help header of "filename"
    ind = find(C2(i,:)==1);
    if any(ind), 
      [cref{3},ind0] = lexsort(M(ind,:)); 
       htmlind{3} = drmask(ind(ind0));
    end %
    
    % Find functions where "filename" is mentioned in help header
    ind=find(C2(:,i)==1);
    if any(ind), 
      [cref{4},ind0] = lexsort(M(ind,:)); 
      htmlind{4} = drmask(ind(ind0));
    end %
    
    
    OC = freadtxt(fullfilename);   % read Original Code
    
    % strings and comments mask  and indices to New lines,...............
    [mask_q,mask_c,inl,linenum] = findquot(OC);
    
    space=' ';
    s2 = OC;
    % Blank the openings of quotes (strings) ........
    nn = find(diff(mask_q)==1)+1;
    s2(nn) = space(ones(size(nn)));
    
    % Blank the comments
    indc = find(mask_c);
    s2(indc) = space(ones(size(indc)));
    s2 = char(s2);

    if 0,
      % Blank the quotes
      string_q = s2;
      nn=find(mask_q);
      string_q(nn) = space(ones(size(nn)));
    end
    %%%%%%%%%% Here we use the TYP control vector %%%%%%%%%%%%
    
    
    if any(typ(1:2)==1),
      %nn2 = find(diff(linenum(indc))>1);
      nn2 = find(diff(indc)>1);
      if isempty(nn2),
	indHH = indc;
      else
	indHH = indc(1:nn2(1)+1); 
      end
      HH = OC(indHH);
      if any(indHH)
	%Blank the leading character of every line in comments
	cMaskHH = mask_c(indHH);
	inlHH = inl(linenum(indHH(1)):linenum(indHH(end)))-indc(1)+1;
	HH([1,inlHH+1]) = space;
      end
    %alternatively
    %  HH = help(fullfilename);
    end
    
    if typ(1)==1, % Find Function Call syntax
      % Should also check on the 1'st line of code
      % contain the function identifier
      [rHH,cHH] = size(HH);
      if rHH >0
	[names1,P0,P1]=getnames(HH);
	res = strmatch('CALL',names1,'exact'); % WAFO style
	
	if ~isempty(res),
	  FCstart = P1(res(1))+2;
	  HHln     = linenum(indHH(FCstart)); % Get the linenumber
	  FCend   = find(linenum(indHH)==HHln);
	  FC  = HH(FCstart:FCend(end));
	  CA{1} = FC;
	end;
      end;
    end;
  
    if typ(2) == 1, % Help Header included
      CA{2} = help(fullfilename);
      if (length(CA{2})< (length(HH)>0))
	CA{2} = HH;
      end
    end
    
    % extract all legal names
    idx=find(typ(3:5));
    if any(idx), 
      % Find any of the following:
      % load-files, print-files or save-files (iz=1,2,3 respectively)
      names2 = getnames(s2);
      for iz=idx
	CA{iz+2} = filesr(names2,iz);
      end
    end
    if typ(6) == 1,% include original m-code
      CA{6} = OC;
    end;

    %%%% The function call to a MAKER file which specifies  %%%%
    %%%% the format used to write out the information about %%%%
    %%%% the m-file.                                        %%%%
    %maker1(htmlind,htmldir,htmlsdir,filename,cref,CA);
    %feval(makername,htmldir,currdir,filename,cref,CA);
    
    [FC,HH,FL,FP,FS,OC]     = deal(CA{:});
    [ref1,ref2,href1,href2] = deal(cref{:});
    tmp = strvcat(ref1,href1);
    [ref1h, ind] = unique(tmp,'rows');
    tmp = cat(2,htmlind{[1 3]});
    ref1hind = tmp(ind);

    %%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%
    fil      = strtok(filename,'.');   % strip the .m extension
    pfile    = fullfile(htmlsdir{ix},[fil '.htm']);
    [fid,msg]=fopen(fullfile(htmldir,pfile),'w');
    if fid==-1,error(msg),end

    fprintf(fid,'<TITLE> %s </TITLE>\n', pfile);
    fprintf(fid,'<BODY>\n');;
    fprintf(fid,'<A HREF = \"index.htm\">[Index for %s]</A>\n',mdir{ix});
    fprintf(fid,'<A HREF = \"../index.htm\">[Return to Master Index]</A>\n');
    fprintf(fid,'<H1>%s</H1>\n' ,filename );
    fprintf(fid,'<H2>(%s)</H2>\n',fullfile(mdir{ix},filename));

    
    isCfile=strcmpi(fil,'contents')|strcmpi(fil,'readme');

    CC = zeros(11,3); % this is the matrix of colors for each
    % call to fntwrt in maker1.m. CC(1,:) is call no.1 and so forth
    % insert the proper color code here.
    
    
    
    % Checking the drive (C:,D:,E: etc):
    dummy = pwd;
    driver = dummy(1:2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%% Writing a header %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    headwrt(fid,[],'Cross-linked m-file',3,0);
%    fntwrt(fid,filename,'CENTER',2,1,CC(1,:));
%    fntwrt(fid,'Located in:','CENTER',1,1,CC(1,:));
%    fntwrt(fid,currdir,'CENTER',1,1,CC(2,:));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fprintf(fid,'%s \n','<P><HR WIDTH="100%"></P>');
    
    %%%%%%%%%%% Writing out the function call %%%%%%%%%%%%%%%%
    if length(FC) >=1 & ~isCfile,
      fntwrt(fid,'Function synopsis','LEFT',2,1,CC(3,:));
      %fprintf(fid,'%s \n','<P><HR WIDTH="100%"></P>');
      fntwrt(fid,FC,'LEFT',1,0,CC(4,:));
    end;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fprintf(fid,'%s \n','<P><HR WIDTH="100%"></P>');
    
    
    %%%% Here we write out the Help Header in the program %%%%%
    %---------------------------------------------------------
    if length(HH) >=1,
      fntwrt(fid,'Help text','LEFT',2,1,CC(5,:));
      totxrefwrt(fid,HH,href1,htmlind{3},htmlsdir,'.htm')
 %     totxrefwrt(fid,HH,currdir,[],ref1h,'.htm');
      %txtwrt(fid,HH);
      fprintf(fid,'%s \n','<P><HR WIDTH="100%"></P>');
    end;
    
    if ~isCfile, % This is a Contents.m or  Readme file -> Stop here
       
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if any(ref1(:))|any(ref2(:))|any(href2(:))
	fprintf(fid,'<HR><H3>Cross-Reference Information</H3>');
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if length(ref1) >=1,  
	fntwrt(fid,['m-files called by ',filename],'LEFT',1,1,CC(6,:))
	
	%%%%% Making a list of pointers to files called 
	[n1,m1]=size(ref1);
	lstwrt(fid,0,0);
	for im =1:n1
	  r2 = strtok(ref1(im,:),'.'); % remove the .m extension
	  % Must add .htm instead of .m since we are
	  % referring to another 
	  r1 = [ r2,'.htm'];
	  nrt = htmlind{1}(im);
	  if nrt~=ix, r1 = fullfile('..',htmlsdir{nrt} ,  r1);end
	  refwrt(fid,r2,r1,1)
	end;
	lstwrt(fid,1,0);
	fprintf(fid,'%s \n','<P><HR WIDTH="100%"></P>');
      end;
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if length(ref2) >=1,
	fntwrt(fid,['m-files that call ',filename],'LEFT',1,1,CC(7,:))
	%%%%% Making a list of pointers to files dependent on current file
	[n2,m2]=size(ref2);
	%r1 = cellstr(ref2);
	%ref1 = strrep(r1,'.m','.htm');
	%ref2 = strrep(r1,'.m','');
	
	lstwrt(fid,0,0);
	for im =1:n2
	  r2 = strtok(ref2(im,:),'.');
	  % Must add .htm instead of .m since we are
	  % referring to another HTML file 
	  r1 = [r2,'.htm'];
	  nrt = htmlind{2}(im);
	  if nrt~=ix, r1 = fullfile('..',htmlsdir{nrt} ,  r1);end
	  refwrt(fid,r2,r1,1)
	end;
	lstwrt(fid,1,0);
	fprintf(fid,'%s \n','<P><HR WIDTH="100%"></P>');
      end;
      
      if length(href2) >=1,
	fntwrt(fid,['m-files that mention ',filename],'LEFT',1,1,CC(7,:))
	%%%%% Making a list of pointers to files where current file is mentioned
	[n2,m2]=size(href2);
	
	lstwrt(fid,0,0);
	for im =1:n2
	  r2 = strtok(href2(im,:),'.'); % remove the .m extension
	  % Must add .htm instead of .m since we are
	  % referring to another HTML file 
	  r1 = [r2,'.htm'];
	  nrt = htmlind{4}(im);
	  if nrt~=ix, r1 = fullfile('..',htmlsdir{nrt},  r1);end
	  refwrt(fid,r2,r1,1)
	end;
	lstwrt(fid,1,0);
	fprintf(fid,'%s \n','<P><HR WIDTH="100%"></P>');
      end;
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      
      %%%% Writing out references to LOADED files      %%%%%%%%
      if length(FL) >=1,
	fntwrt(fid,['Files that are LOADED into ',filename],'LEFT',1,1,CC(8,:))
	[n,m]=size(FL);
	lstwrt(fid,0,0);
	for im =1:n
	  inpput_file_names =[FL(im,:) blanks(1)];
	  tmp = deblank(inpput_file_names);
	  tmp = strrep(tmp,driver,''); %supposed to remove the driver
	  aa = findstr(tmp,'.');
	  if length(aa) <1,
	    tmp = [tmp,'.mat'];
	  end;
	  if tmp(1) == '\'
	    tmp2 = ['file:///',driver(1),'|' tmp];
	  else
	    tmp2 = ['file:///',driver(1),'|' '/' tmp];
	  end;
	  tmp2 = strrep(tmp2,'\','/');
	  refwrt(fid,tmp,tmp2,1);
	end;
	lstwrt(fid,1,0);
	fprintf(fid,'%s \n','<P><HR WIDTH="100%"></P>');
      end;
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      
      %%%% Writing out references to PRINTED files     %%%%%%%%
      if length(FP) >=1,
	fntwrt(fid,['Files that are PRINTED from ',filename],'LEFT',1,1,CC(9,:))
	[n,m]=size(FP);
	lstwrt(fid,0,0);
	for im =1:n
	  printed_file_names = [FP(im,:) blanks(1)];
	  tmp = deblank(printed_file_names);
	  tmp = strrep(tmp,driver,''); %supposed to remove the driver
	  aa = findstr(tmp,'.');
	  if length(aa) <1,
	    tmp = [tmp,'.mat'];
	  end;
	  if tmp(1) == '\'
	    tmp2 = ['file:///',driver(1),'|' tmp];
	  else
	    tmp2 = ['file:///',driver(1),'|' '/' tmp];
	  end;
	  tmp2 = strrep(tmp2,'\','/');
	  refwrt(fid,tmp,tmp2,1);
	end;
	lstwrt(fid,1,0);
	fprintf(fid,'%s \n','<P><HR WIDTH="100%"></P>');
      end;
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

      %%%% Writing out references to SAVED files       %%%%%%%%
      if length(FS) >=1,
	fntwrt(fid,['Files that are SAVED from ',filename],'LEFT',1,1,CC(10,:))
	[n,m]=size(FS);
	lstwrt(fid,0,0);
	for im =1:n
	  saved_file_names = [currdir '/' FS(im,:) blanks(1)];
	  tmp = deblank(saved_file_names);
	  tmp = strrep(tmp,driver,''); %supposed to remove the driver
	  aa = findstr(tmp,'.');
	  if length(aa) <1,
	    tmp = [tmp,'.mat'];
	  end;
	  if tmp(1) == '\'
	    tmp2 = ['file:///',driver(1),'|' tmp];
	  else
	    tmp2 = ['file:///',driver(1),'|' '/' tmp];
	  end;
	  tmp2 = strrep(tmp2,'\','/');
	  refwrt(fid,tmp,tmp2,1);
	end;
	lstwrt(fid,1,0);
	fprintf(fid,'%s \n','<P><HR WIDTH="100%"></P>');
      end;
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if length(OC) >=1,
	fntwrt(fid,['All the cross references in the source m-code of ',filename],'LEFT',1,1,CC(11,:))
	%totxref(fid,currdir,[],filename,ref_to,'.htm');
	totxrefwrt(fid,OC,ref1h,ref1hind,htmlsdir,'.htm')
%	totxrefwrt(fid,OC,currdir,[],ref1h,'.htm');
	fprintf(fid,'%s \n','<P><HR WIDTH="100%"></P>');
      end;
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    sigwrt(fid);
    fclose(fid);
  
  end  
end;  
return

function C = filesr(names,typ)
% C = filesr(names,typ)
% Returns the names of files read or written by
% a selected m-file
%
% C      = a string matrix where each row contains the name of the
%          specific file type. Extensions will be included as they
%          are specified in the loading/saving/printing commands.
% names  = character array with all the legal names of the m-file 
%          to be searched. 
% type   = 1 : load-files   (look for the "load" filename)
%          2 : print-files (look for print -dxxxx filename)
%          3 : save-files  (look for "save" filename)
%          4 : fopen-files (look for "fopen" filename)
%
% OBS! Only one command per line in m-file!
%

error(nargchk(2,2,nargin))
C =[]; 
%names = getnames(names);
nr = size(names,1); % number of rows
if nr<1,return; end

switch typ ,
  case 1, %%%%%%%%%%%% Look for LOAD   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ix = strmatch('load',names,'exact');
    if any(ix),
      ix(ix>=nr)=[];
      C=names(ix+1,:); % the word following load is the filename!
    end
  case 2, %%%%%%%%%%%% Look for print   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ix = strmatch('print',names,'exact');
    if ~any(ix), return, end
    ix(ix>=nr-1)=[];
    dev = names(ix+1); % device options
    ind=find(isprntoption(dev));
    if ~any(ind), 
      disp(['1)Can not handle print options  in filsr' num2str(length(ind))]),
    end
    C1=names(ix+2,:); % the word following the print device is the filename!
    ind=find(isprntoption(C1));
    if ~any(ind),
      disp(['2) Can not handle print options in filsr' num2str(length(ind))]),
    end
    C2=cellstr(C1);
    for ix=1:length(C2)
      ext = prntext(dev(ix,:));
      if isempty(ext),
	disp('Invalid print device')
      else
	C2(ix)=  {[C2{ix} ext]};
      end
    end
    C=char(C2);
    
  case 3,%%%%%%%%%%%% Look for save   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ix = strmatch('save',names,'exact');
    if any(ix)
      ix(ix>=nr)=[];
      C=names(ix+1,:); % the word following save is the filename!
    end
  case 4,%%%%%%%%%%%% Look for fopen   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ix = strmatch('fopen',names,'exact');
    if any(ix)
      ix(ix>=nr)=[];
      C=names(ix+1,:); % the word following save is the filename!
    end
end; %switch typ
return



function c  = prntext(pr0)
% PRNTEXT Return correct file extension given the print device  
%
% CALL: [c,flag] = prntext(device);
%
% device  = 'deps', 'deps2','depsc', 'depsc2', 'dps', 'dps2','dpsc', 
%           'dpsc2' or 'dhpgl'   
% c       = string containing the correct extension,
%           i.e. '.eps','.ps', or '.hgl'                  
%           (encapsulated postscript, postscript or HP plotter,
%               respectively)

ss = deblank(pr0);
switch ss
  case {'deps','depsc','deps2','depsc2'}, c = '.eps'; % Encapsulated PostScript
  case {'dps','dpsc','dps2','dpsc2'},     c = '.ps' ; % PostScript
  case 'dhpgl',                           c = '.hgl';
  case 'ill',                             c = '.ai' ;
  case 'mfile',                           c = '';
  case 'jpeg',                            c = '.jpg';
  case {'tiff',  'tiffnocompression'},    c = '.tif';
  case  'meta',                           c = '.emf';
  case 'bitmap',                          c = '.bmp';
  case 	{'laserjet','ljetplus','ljet2p','ljet3','cdeskjet','cdjcolor','cdjmono','deskjet','cdj550','djet500','paintjet','pjetxl','pjxl','pjxl300','cdj500','dnj650c','bj10e','bj200','bjc600'}, c = '.jet';
  case {'epson','epsonc','eps9high'},   c = '.ep';  
  case 'ibmpro',  c = '.ibm';
  case 'ln03',  c = '.ln3';
  case {'pcxmono','pcx16','pcx256','pcx24b'}, c ='.pcx';
  case {'bmp16m','bmp256'},      c= '.bmp';
  case {'pbm','pbmraw'},       c = '.pbm';
  case {'pgm','pgmraw'},       c = '.pgm';
  case {'ppm','ppmraw'},       c = '.ppm';
  case 'sgirgb',       c = '.sgi';
    otherwise, c=[]; % unknown device
end;
return
	
       
function flag=isprntoption(str)
% ISPRNTOPTION  returns ones for the words that are legal print options
%
%  CALL: flag = isprntoption(words)
%

% legal options to the printer
options = strvcat('v', 'loose','epsi', 'tiff', 'append','adobecset','cmyk',...
    'r', 'noui', 'painters', 'zbuffer', 'DEBUG');
[n m ]=size(str);
flag=zeros(n,1);

for ix=1:n,
  flag(ix) = any(strmatch(deblank(str(ix,:)),options,'exact'));
end
return




function [] = headwrt(fid,directory,overhead,fontsize,type)
% HEADWRT  Write a simple header for a HTML document to file.
%
% CALL: headwrt(fid,directory,overhead,fontsize,type)
%
% fid       = pointer to current HTML-file
% directory = current directory
% overhead  = string with the title of the HTML doc.
% fontsize  = size of title string
% type      = 0-> bold, 
%             1-> italics
% 
%  Note: uses the global variable HTINI

% History:
% revised pab 17.10.2000
%  - changed name HEAD_mak -> headwrt
%  - added global HTINI instead of iniread1-3
% Copyright (c) 1996 B.K. Alsberg
%
error(nargchk(3,5,nargin))
if nargin<4|isempty(fontsize), fontsize = 3; end
if nargin<5|isempty(type),    type = 1; end

global HTINI
gif_file = HTINI.bg;       % background
dirloc   = HTINI.gif_home; % location of GIF's
color    = HTINI.color_idx;

%%%%%%%% First we write the HTML header info %%%%%%%%%
s(1,:) = '<HTML>                                             ';
s(2,:) = '<HEAD>                                             ';
s(3,:) = '  <TITLE>/</TITLE>                                 ';
s(4,:) = '  <META NAME="GENERATOR" CONTENT="HTML-TOOLBOX ">  ';
s(5,:) = '</HEAD>                                            ';
s(6,:) = '<BODY>                                             ';
s(7,:) = '</BODY>                                            ';
s(8,:) = '</HTML>                                            ';

for j = 1:5
 fprintf(fid,'%s \n',s(j,:));
end;

%% Here we print out the code of selecting a background GIF-file:
if length(gif_file) >0,
 background = ['<BODY BACKGROUND="',dirloc,gif_file,'">'];
 fprintf(fid,'%s \n',background); 
end;

fprintf(fid,'%s \n','<HR WIDTH="100%"></P>');

fntwrt(fid,overhead,'CENTER',fontsize,type,color);
fntwrt(fid,directory,'CENTER',fontsize-2,type,color);
return

function  refwrt(fid,refstring,filename,lst,target)
% REFWRT  Writes a HTML reference to a specific file name
%
% CALL: refwrt(fid,refstring,filename,lst,target);
%
% fid               = file pointer to the file you
%                     are currently working with
% refstring         = the string in the HTML which
%                     identifies your link (i.e. the
%                     underlined writing)
% filename          = the file your are creating a link to
% lst               = 0-> not member of a bullet list
%                     1-> member of a bullet list
% target            = use only for frames. Indicate the TARGET
%                     frame where the html file will be displayed



%History
% revised pab 17.10.2000
%

%if lst == 1,
% str = ['<LI><A HREF ="',filename,'">',refstring,'</A></LI>']; 
% fprintf(fid,'%s \n',str);
%elseif lst == 0
% str = ['<A HREF ="',filename,'">',refstring,'</A>']; 
% fprintf(fid,'%s \n',str);
%end;

if nargin == 5,
 trg = [' TARGET = ' '"' target '"'];
else
 trg = [];
end;

if lst == 1,
  str = ['<LI><A HREF ="',filename,'"' trg  '>',refstring,'</A></LI>']; 
else %if lst == 0
  str = ['<A HREF ="',filename,'"'  trg  '>',refstring,'</A>']; 
end;
fprintf(fid,'%s \n',str);
return


function sigwrt(fid)
% SIGWRT Write signature to HTML-file
% 
% CALL: sigwrt(fid)
%
% fid      = is an integer file identifier obtained from FOPEN.
% 
%  SIGWRT Reads the global variable HTINI with fields:
% Name     = string with user name
% par      = parameter vector:
%  par(1)  = 0-> do not include the date, 1-> include
%  par(2)  = 0-> do not include the time, 1-> include
%  par(3)  = 0-> do not include the copyright GIF, 1-> include
% email    = string containing the e-mail address. Optional
%

% History:
% revised pab 17.10.2000
%  changed name from sig_wrt1 > sigwrt
% Copyright (c) 1996 B.K. Alsberg
%


global HTINI
%% Initialize the parameters to the signature file by %%%%%
%% reading from the global variable HTINI
Name  = HTINI.name;
email = HTINI.email;
par   = HTINI.par;


%%%% Init of color,font,size for first and second part of signature %%%
color_sig1  = HTINI.color_sig1;
color_email = HTINI.color_email;
font        = HTINI.font;
fontsize    = HTINI.fontsize;
gif_home    = HTINI.gif_home;
cp          = HTINI.cp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if par(1) == 1
  date = datestr(now,1);%[day,'/',month,' ',year];
else
  date = [];
end;

if par(2) == 1,
  time = datestr(now,15);%[hour,':',minute];
else
  time =[];
end;


if par(3) == 1,
  image_str = ['<IMG SRC = "',fullfile('..',cp),'" >'];
else
  image_str =[];
end;


if length(email) >0,
  email_str = ['E-mail:</B><A HREF = "mailto:',email,'">',email,'</A></H4>'];
else
  email_str =  [];
end;

message = 'Written by ';
str1    = [message,' ',Name,' ',time,' ',date,' ',image_str];

%here we transform from [0,1] space into [0,255] space
color_sig2   = round(lintrans(color_sig1,0,1,0,255));
color_email2 = round(lintrans(color_email,0,1,0,255));


if isempty(font)|font == 0
  front = '<B>';
  back  = '</B>';
elseif font ==1,
  front = '<I>';
  back  = '</I>';
else
  front = '';
  back  = '';
end

if fontsize >0
  fnt = ['+',int2str(fontsize)];
else
  fnt = ['-',int2str(fontsize)];
end;

fz = '<FONT SIZE=';
p1 = '<P>';
p2 = '</P>';
fnt2 = '</FONT>';
h0 = '>';

%%%%%% Here we convert the color vector into proper HTML code %%%
red = dec2hex(color_sig2(1));
if length(red) == 1, red = ['0',red]; end;
green = dec2hex(color_sig2(2));
if length(green) == 1, green = ['0',green]; end;
blue = dec2hex(color_sig2(3));
if length(blue) == 1, blue = ['0',blue]; end;
htcd1 = ['#',red,green,blue];
fntcol1 = ['<FONT COLOR="',htcd1,'">'];
clear red green blue
%%%%%%%%%%%%%%%%%%%%%%%%%%% Second line %%%%%%%%%%%%%%%%%%%

red = dec2hex(color_email2(1));
if length(red) == 1, red = ['0',red]; end;
green = dec2hex(color_email2(2));
if length(green) == 1, green = ['0',green]; end;
blue = dec2hex(color_email2(3));
if length(blue) == 1, blue = ['0',blue]; end;
htcd2 = ['#',red,green,blue];
fntcol2 = ['<FONT COLOR="',htcd2,'">'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                                                                                                   
output1 = [p1,front,fntcol1,fz,fnt,h0,str1,fnt2,fnt2,back,p2];
fprintf(fid,'%s \n',output1);

output2 = [p1,front,fntcol2,fz,fnt,h0,email_str,fnt2,fnt2,back,p2];
fprintf(fid,'%s \n',output2);
return

function [ut] = lintrans(x,x1,x2,y1,y2)
% LINTRANS  Linear transformation of variables
%
% CALL: [y] = lintrans(x,x1,x2,y1,y2)
%
%  y  = transformed values
%  x  = original values
%
% (x1,x2) ---> (y1,y2). 
% Example: Map x = 1:4 into y = -5:5 :
%
% y = lintrans(x,1,4,-5,5);


a = (y1-y2)/(x1-x2);
b = y1 - a*x1;
ut = a*x + b;
return


function txtwrt(fid,U)
% TXTWRT  Writes the text in U in file fid 
%
% CALL: txtwrt(fid,U)
%

[n,m]=size(U);

fprintf(fid,'%s \n','<pre>');
for i =1:n
 fprintf(fid,'%s \n',deblank(U(i,:)));
end;
fprintf(fid,'%s \n','</pre>');

return
function totxrefwrt(fid,str,refm,ind,htmlsdir,ext)
% TOTXREFWRT Write string with inserted links to other m-files to HTML-file
%
% CALL: totxrefwrt(fid,str,refm,ind,htmlsdir,ext)
%
% fid          = file pointer to the current HTML file
% str          = string to be searced for crossreferences.
% refm         = string matrix of all functions to cross reference.
% ind          = indices to html subdirectories.
% htmlsdir     = cellarray of html subdirectories where refm files are located.
% ext          = extension for the files in refm (default .htm)

% History
% -revised pab 15.05.2001
% - revised pab 17.10.2000
% changed name from totxref -> totxrefwrt
% Copyright (c) 1996 B.K. Alsberg
%

error(nargchk(5,6,nargin))
if nargin<6|isempty(ext), ext ='.htm';end
%if nargin<4|isempty(htmlsdir)), htmlsdir =[];end
[n2,m2]=size(refm);

% extract all legal names from string
[names,P0,P1] = getnames(str); 
names2 = lower(names);

href = [];no = [];nc=[];
for j = 1:n2,
  %%%% Here we remove the extension     %%%%%%%%%
  refstring = deblank(strtok(refm(j,:),'.'));
  iy = strmatch(refstring,names2,'exact');
  if any(iy),
    Niy = length(iy);
    
    filepointer = [refstring,ext];
    nrt = ind(j);
    if ~isnan(nrt), 
      filepointer = fullfile('..',htmlsdir{nrt},filepointer);
    end;
    %%% Here we construct the HTML string used in file links: %%
    href0 = ['<A HREF ="',filepointer,'">'];%,refstring,'</A>'];
    for iz =iy(:).'
      href = strvcat(href,[href0,deblank(names(iz,:)),'</A>']);
    end
    no = cat(2,no,P0(iy));
    nc = cat(2,nc,P1(iy));
  end
end

str = insert(str,href,no,nc+1);
 
fprintf(fid,'%s \n','<pre>');
fprintf(fid,'%s \n',str);
fprintf(fid,'%s \n','</pre>');

return

function maker1(htmldir,currdir,filename,cref,CA)
% MAKER1 Generates the crosslinked HTML files
%
% CALL: maker1(htmldir,currdir,filename,cref,CA)
%
% htmldir  = root directory for the html files
% currdir  = directory containing the m-file
% filename = name of m-file (without directory) we want to 
%            describe in a HTML file. Must include the .m extension
% cref     = {ref1, ref2, href1, href2} cell array of string matrices
%            with crossreference information:
%              ref1  : m-files called by "filename" 
%              ref2  : m-files dependent on "filename" 
%              href1 : m-files mentioned in the header of "filename" 
%              href2 : m-files where "filename" is mentioned in the header
% CA       = {FC, HH, FL, FP, FS, OC} a cell array of character arrays 
%            to include in the HTML file (size [1 6]) where:
%               FC : Function Call syntax of the m-file 
%               HH : Help Header
%               FL : Files Loaded into m-file
%               FP : Files Printed from m-file
%               FS : Files Saved from m-file
%               OC : Original m-Code.
%            Note: Empty strings will not be included in the HTML-file.


% History
% revised pab 18.10.2000
% Copyright (c) 1996 B.K. Alsberg
%
error(nargchk(5,5,nargin))

[FC,HH,FL,FP,FS,OC]     = deal(CA{:});
[ref1,ref2,href1,href2] = deal(cref{:});
tmp = strvcat(ref1,href1);
ref1h = munique(tmp);


%%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%
fil = strtok(filename,'.');   % strip the .m extension
fid = fopen([htmldir,filesep,fil,'.htm'],'w');
%fid = fopen([fil,'.htm'],'w');

isCfile=strcmpi(fil,'contents');
isCfile=strcmpi(fil,'readme')|isCfile;

CC = zeros(11,3); % this is the matrix of colors for each
% call to fntwrt in maker1.m. CC(1,:) is call no.1 and so forth
% insert the proper color code here.



% Checking the drive (C:,D:,E: etc):
dummy = pwd;
driver = dummy(1:2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%% Writing a header %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 headwrt(fid,[],'Cross-linked m-file',3,0);
 fntwrt(fid,filename,'CENTER',2,1,CC(1,:));
 fntwrt(fid,'Located in:','CENTER',1,1,CC(1,:));
 fntwrt(fid,currdir,'CENTER',1,1,CC(2,:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(fid,'%s \n','<P><HR WIDTH="100%"></P>');

%%%%%%%%%%% Writing out the function call %%%%%%%%%%%%%%%%
if length(FC) >=1 & ~isCfile,
 fntwrt(fid,'Function synopsis','LEFT',2,1,CC(3,:));
 %fprintf(fid,'%s \n','<P><HR WIDTH="100%"></P>');
 fntwrt(fid,FC,'LEFT',1,0,CC(4,:));
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(fid,'%s \n','<P><HR WIDTH="100%"></P>');


%%%% Here we write out the Help Header in the program %%%%%
%---------------------------------------------------------
if length(HH) >=1,
 fntwrt(fid,'Function comments','LEFT',2,1,CC(5,:));
 totxrefwrt(fid,HH,currdir,[],ref1h,'.htm');
 %txtwrt(fid,HH);
 fprintf(fid,'%s \n','<P><HR WIDTH="100%"></P>');
end;

if isCfile, % This is a Contents.m or  Readme file -> Stop here
  sigwrt(fid);
  fclose(fid);
  return
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(ref1) >=1,  
 fntwrt(fid,['m-files called by ',filename],'LEFT',1,1,CC(6,:))
 %fprintf(fid,'%s \n','<P><HR WIDTH="100%"></P>');
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %%%%% Making a list of pointers to files called 
 [n1,m1]=size(ref1);
 lstwrt(fid,0,0);
 for i =1:n1
  %r1 = deblank(ref1(i,:));
  %r2 = r1(1:length(r1)-2);
  r2 = strtok(ref1(i,:),'.'); % remove the .m extension
  % Must add .htm instead of .m since we are
  % referring to another 
  r1 = [r2,'.htm'];
  refwrt(fid,r2,r1,1)
 end;
 lstwrt(fid,1,0);
 fprintf(fid,'%s \n','<P><HR WIDTH="100%"></P>');
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(ref2) >=1,
 fntwrt(fid,['m-files that call ',filename],'LEFT',1,1,CC(7,:))
 %%%%% Making a list of pointers to files dependent on current file
 [n2,m2]=size(ref2);
 %r1 = cellstr(ref2);
 %ref1 = strrep(r1,'.m','.htm');
 %ref2 = strrep(r1,'.m','');
 
 lstwrt(fid,0,0);
 for i =1:n2
%  r1 = deblank(ref2(i,:));
%  r2 = r1(1:length(r1)-2);
  r2 = strtok(ref2(i,:),'.');
  % Must add .htm instead of .m since we are
  % referring to another HTML file 
  r1 = [r2,'.htm'];
  refwrt(fid,r2,r1,1)
 end;
 lstwrt(fid,1,0);
 fprintf(fid,'%s \n','<P><HR WIDTH="100%"></P>');
end;

if length(href2) >=1,
 fntwrt(fid,['m-files that mention ',filename],'LEFT',1,1,CC(7,:))
 %%%%% Making a list of pointers to files where current file is mentioned
 [n2,m2]=size(href2);
 %r1 = cellstr(ref2);
 %ref1 = strrep(r1,'.m','.htm');
 %ref2 = strrep(r1,'.m','');
 
 lstwrt(fid,0,0);
 for i =1:n2
%  r1 = deblank(href2(i,:));
%  r2 = r1(1:length(r1)-2);
  r2 = strtok(href2(i,:),'.'); % remove the .m extension
  % Must add .htm instead of .m since we are
  % referring to another HTML file 
  r1 = [r2,'.htm'];
  refwrt(fid,r2,r1,1)
 end;
 lstwrt(fid,1,0);
 fprintf(fid,'%s \n','<P><HR WIDTH="100%"></P>');
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% Writing out references to LOADED files      %%%%%%%%
if length(FL) >=1,
 fntwrt(fid,['Files that are LOADED into ',filename],'LEFT',1,1,CC(8,:))
 [n,m]=size(FL);
 lstwrt(fid,0,0);
 for i =1:n
   inpput_file_names =[FL(i,:) blanks(1)];
   tmp = deblank(inpput_file_names);
   tmp = strrep(tmp,driver,''); %supposed to remove the driver
   aa = findstr(tmp,'.');
   if length(aa) <1,
    tmp = [tmp,'.mat'];
   end;
   if tmp(1) == '\'
    tmp2 = ['file:///',driver(1),'|' tmp];
   else
    tmp2 = ['file:///',driver(1),'|' '/' tmp];
   end;
   tmp2 = strrep(tmp2,'\','/');
   refwrt(fid,tmp,tmp2,1);
 end;
 lstwrt(fid,1,0);
 fprintf(fid,'%s \n','<P><HR WIDTH="100%"></P>');
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% Writing out references to PRINTED files     %%%%%%%%
if length(FP) >=1,
 fntwrt(fid,['Files that are PRINTED from ',filename],'LEFT',1,1,CC(9,:))
 [n,m]=size(FP);
 lstwrt(fid,0,0);
 for i =1:n
   printed_file_names = [FP(i,:) blanks(1)];
   tmp = deblank(printed_file_names);
   tmp = strrep(tmp,driver,''); %supposed to remove the driver
   aa = findstr(tmp,'.');
   if length(aa) <1,
    tmp = [tmp,'.mat'];
   end;
   if tmp(1) == '\'
    tmp2 = ['file:///',driver(1),'|' tmp];
   else
    tmp2 = ['file:///',driver(1),'|' '/' tmp];
   end;
   tmp2 = strrep(tmp2,'\','/');
   refwrt(fid,tmp,tmp2,1);
 end;
 lstwrt(fid,1,0);
 fprintf(fid,'%s \n','<P><HR WIDTH="100%"></P>');
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% Writing out references to SAVED files       %%%%%%%%
if length(FS) >=1,
 fntwrt(fid,['Files that are SAVED from ',filename],'LEFT',1,1,CC(10,:))
 [n,m]=size(FS);
 lstwrt(fid,0,0);
 for i =1:n
   saved_file_names = [currdir '/' FS(i,:) blanks(1)];
   tmp = deblank(saved_file_names);
   tmp = strrep(tmp,driver,''); %supposed to remove the driver
   aa = findstr(tmp,'.');
   if length(aa) <1,
    tmp = [tmp,'.mat'];
   end;
   if tmp(1) == '\'
    tmp2 = ['file:///',driver(1),'|' tmp];
   else
    tmp2 = ['file:///',driver(1),'|' '/' tmp];
   end;
   tmp2 = strrep(tmp2,'\','/');
   refwrt(fid,tmp,tmp2,1);
 end;
 lstwrt(fid,1,0);
 fprintf(fid,'%s \n','<P><HR WIDTH="100%"></P>');
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(OC) >=1,
 fntwrt(fid,['All the cross references in the source m-code of ',filename],'LEFT',1,1,CC(11,:))
 %totxref(fid,currdir,[],filename,ref_to,'.htm');
 totxrefwrt(fid,OC,currdir,[],ref1h,'.htm');
 fprintf(fid,'%s \n','<P><HR WIDTH="100%"></P>');
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sigwrt(fid);
fclose(fid);
