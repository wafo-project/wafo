#!/usr/local/bin/perl

$URL = "http://www.nd.edu/Courses/kantor/matlab/mat2html";

@keyWords = ('break', 'else', 'elseif', 'end', 'for', 'if', 'return',
				 'while', 'function', 'global', 'eval', 'feval', 'nargchk', 'pause', 'menu',
				 'keyboard', 'input', 'ans', 'computer', 'eps', 'flops', 'inf', 'NaN', 
				 'nargin', 'nargout', 'pi', 'realmax', 'realmin');

# READ THE HELP SECTION BELOW FOR INFORMATION ON HOW TO USE THIS SCRIPT.
#
# To install, simply save this function at mat2html, and make it 
# executible with
#
#     chmod +x mat2html
#
# This requires a perl interpreter and access to the standard perl
# library. Do not contact me about installing perl, I will ignore
# you. Instead, please check out the following standard sites:
#
#    ftp.uu.net                      gnu/perl*
#    jpl-decvax.jpl.nasa.gov         pub/perl*
#    archive.cis.ohio-state.edu      perl/
#
#
#                13 Jan 1995: Fixed a bug in tracking cross-
#                             references, and other refinements.
# Version 1.0 -> 12 Jan 1995: Major workover to allow index of
#                             multiple matlab directories.
#                10 Jan 1995: Added cross-links to the help text
#                 9 Jan 1995: Looks for a Readme.m file to put
#                             to put in the directory index.
#                 9 Jan 1995: Fixed html cross-links
#
# First posted to comp.soft-sys.matlab on January 9, 1995

# Jeffrey C. Kantor
# kantor.1@nd.edu

# Copyright (c) Jeffrey C. Kantor 1995.
# All rights reserved
#
# Please feel free to use this script under the conditions of the
# standard GNU public license.

######################################################################
#
#  Subroutines
#
######################################################################

# From the standard Perl Library

require 'getopts.pl';
require 'ctime.pl';

# Translate special html characters

sub htmlchar {
  s/&/&amp;/g;
  s/\</&lt;/g;
  s/\>/&gt;/g;
  s/\"/&quot;/g;
}

# Write an html anchor, and keep track of how many have been
# written using the global variable $anchor_count

$anchor_count = 0;

sub anchor_url{
  local($url,$name) = @_;
  $anchor_count++;
  return("<A HREF = \"$url\">$name<\/A>");
}

# A standard tagline is added to each page created.
# The total page count is kept in $page_count

$page_count = 0;

sub tagline {
  $page_count++;
  print HFILE "Produced by ",&anchor_url($URL,'mat2html')," on $date<BR>\n";
  print HFILE "Cross-Directory links are: ";
  if ($opt_g) {print HFILE "ON<BR>\n"} else {print HFILE "OFF<BR>\n"};
}

# Assume MFILE and HFILE are open, $headline contains a section
# title, and @zref contains cross-reference names

sub writehelpblock {
  do {$_ = <MFILE>} until /^\s*%/ || eof;
  if (!eof) {
    print HFILE "<HR><H3>$headline</H3>\n";
    print HFILE "<pre>\n";
    while (/^\s*%/) {
      s/^\s*%//;
      &htmlchar;
      foreach $z (@zref) {
        next if /<A.*$z.*A>/;
        $html = "../$hdir{$mdir{$mfile{$z}}}/$z.html";
        s/(\W+)($z)(\W+)/$1.&anchor_url($html,$2).$3/eig;
      }
      print HFILE $_;
      $_ = <MFILE>;
    }
    print HFILE "</pre>\n";
  }
}

######################################################################
#
#  Command line processing
#
######################################################################

# Get the current date string

$date = &ctime(time);

# Read the command line arguments

&Getopts('ghipqH:M:');

warn("WARNING: Options -M and -p are incompatable.\n") if ($opt_M && $opt_p);

if (($opt_h) || ($opt_M && $opt_p)) {
  print <<EOF;

Usage:

  mat2html [-i] [-q] [-g] [-M matlab_dir] [-H html_dir]
  mat2html [-i] [-q] [-g] [-p] [-H html_dir]

mat2html reads a list of matlab .m files and/or directories from the
standard input to produce a hypertext documentation suitable for
browsing with a WWW browser such as Mosaic or Netscape.  An file
html_dir/index.html in written. Subdirectories are written for matlab
directory that is encountered containing
.html files corresponding to every .m file.

Help:

  -h   Print this help message and exit.

  -q   Quiet operation. Default is verbose.

Output Options:

  -H  Specify a directory to place the html files and subdirectories.
      The default is the current directory. If necessary, an attempt
      is made to create the directory. The file index.html is placed
      in this directory.

  -i  Include matlab source code in the html documentation

  -g  Do global hypertext links, that is, hypertext links among
      separate matlab directories. Default is to do hypertext
      links only among functions in the same directory.

Matlab Source Options:

      The standard input is looked at first. If there is nothing there,
      then we look in the current directory for .m files as if we did a
      -M . option.

  -M  Specify a root matlab directory to search. The standard input
      is not read.

  -p  Search the matlab path obtained by matlab -n. Options -M
      and -p are incompatable.

Typical usages are:

  mat2html 
      Produces a file index.html in the current directory, and a 
      .html file corresponding to each .m file in the current directory.
      An index.html file is produced in the current directory.

  mat2html -M matlab -H html
      Same as above, except it looks in subdirectory matlab of
      the current directory for the .m files, and puts the results
      in subdirectory html. Useful for separating the .m files
      from the documentation.

  mat2html -p -H html
      Creates documentation for everything on the default matlab
      path, and puts it into a directory html.

  ls *.m | mat2html -H html
      Index the .m files in the current directory

  find . -name "*.m" -print | mat2html -H html -i
      The find command recursively builds a list of all .m files
      to be found the directory and its subdirectories. These
      are then processed and the html files put in the directory
      html. The matlab source code is included in the html files.

EOF
  exit;
}

# Read arguments and initialize

if ($opt_H) {$hroot = $opt_H} else {$hroot = "."};
if ($opt_q) {$verbose = 0} else {$verbose = 1};

######################################################################
#
#  Reading the input
#
######################################################################

# Get the list of files and directories to examine, put in @xfiles
# Start by checking the standard input.

undef(@xfiles);

if (!-t) {

  # STDIN is connected to some input
  print "Reading from the standard input\n" if $verbose;

  warn("Reading from STDIN, ignoring redundant -M option.\n") if ($opt_M);

  while (<>) {
    chop;
    @files = split;
    @xfiles = (@xfiles,@files);
  }

} elsif ($opt_p) {

  # Get the matlab path
  print "Running matlab -n to get the matlab path.\n" if $verbose;

  open(MATLAB,"matlab -n|") || die("Can't run matlab -n to get the path");
  do {$_ = <MATLAB>} until /MATLABPATH/ || eof;
  while (<MATLAB>) {
    last if /^----/;
    chop;
    s/\s*//g;
    push(@xfiles,$_);
  }
  close(MATLAB);

} elsif ($opt_M) {

  # We need to find the matlab directories to search
  print "Setting matlab directory to $opt_M\n" if $verbose;

  if (-d $opt_M) {
    @xfiles = ($opt_M);
  } else {
    die("Specified -M $opt_M is not a directory\n");
  }

} else {

  print "Setting matlab directory to the current directory.\n" if $verbose;
  @xfiles = (".");

}

# Now process the list of @xfiles to get a list of .m files @mfiles

print "Finding .m files in the search path.\n" if $verbose;

undef(@mfiles);
foreach (@xfiles) {

  # chop off any trailing separators, sometimes added by ls

  s/\/$//;

  # Add to mfiles if it exits and text file with a .m suffix

  if ((/\.m$/) && (-e) && (-T)) {
    push(@mfiles,$_);
    next;
  }

  # If it's a directory, then read it

  if (-d) {
    opendir(MDIR,$_) || die("Cannot open directory $_\n");
    @files = grep(/\.m$/,readdir(MDIR));
    foreach $file (@files) {
      push(@mfiles,"$_/$file");
    }
  }
}

$n = $#mfiles + 1;
print "Found $n .m files.\n" if $verbose;

undef($n);
undef(@xfiles);
undef($mroot);

######################################################################
#
#  Parse the matlab file names
#
######################################################################

# Now we need to parse the mfile names to obtain for each file
#  $name{$file} - a matlab identifier used to search the texts
#  $mdir{$file} - the directory in which the file is found

print "Parsing the matlab file names.\n" if $verbose;

undef(%mdir);
foreach (@mfiles) {
  local($x) = $_;
  $x =~ s/\.m$//;
  split("/",$x);
  $name{$_} = pop(@_);
  $mdir{$_} = join('/',@_);
  $mdir{$_} = "." if $mdir{$_} eq "";
}

# Compute a unique list of matlab identifier names, put in @names.

grep($count{$_}++,values(%name));
@names = sort(keys(%count));
undef(%count);

$n = $#names + 1;
print "Found $n unique matlab identifiers.\n" if $verbose;
undef($n);

# Now we have a problem. Each matlab name may be associated with more
# than one .m file. The order of @mfiles is the order they would be
# encountered in a standard matlab search, so we preserve that order.
# We compute:
#
# $mfile{$name}  - look up the first mfile associated with a $name
#
# Confused yet?

print "Linking each matlab identifier to a unique .m file.\n" if $verbose;

foreach $n (@names) {
  foreach (@mfiles) {
    $file = $_;
    last if ($n eq $name{$_});
  }
  $mfile{$n} = $file;
}

# Compute the set of unique matlab directory names, put in @mdirs

grep($count{$_}++,values(%mdir));
@mdirs = sort(keys(%count));
undef(%count);

$n = $#mdirs + 1;
print "Found $n unique matlab directories in the search path.\n" if $verbose;
undef($n);

######################################################################
#
#  Read the m-files, and compute cross-reference information
#
######################################################################

# Read each file and tabulate the distinct alphanumeric identifiers in
# an array of symbols. This is used to compute cross references within
# a matlab directory. Also scan for:
#
#   synposis  - The function declaration line
#   lookfor   - The first line of the help text
#   mtype     - File type, either "function" or "script"
#   cref,href - Arrays o cross references
#
# Dependency matrix value $cref{$x,$y} is 1 if $x refers to $y in the
# code section. $href{$x,$y} is 1 if $x refers to $y in the help
# lines.
  
print "Read and compute cross-references among the .m files.\n" if $verbose;

foreach $file (@mfiles) {

  open(MFILE,"<$file") || die("Cannot open $file\n");

  print "Reading $file\n" if $verbose;

  while (<MFILE>) {
    chop;

    # If it's the function declaration line, then store it and skip
    if (/^\s*function/) {
      s/^\s*function\s*//;
      $synopsis{$file} = $_;
      $mtype{$file} = "function";
      next;
    }

    # Compress multiple %'s to a single %
    s/%+/%/g;

    # Process comment lines and code lines separately

    if (/^\s*%/) {

      # cut out comment marker and surrounding white space
      s/^\s*%\s*//;

      # Store first comment line in lookfor
      if (!$lookfor{$file}) {
        $lookfor{$file} = $_;
      }

      # Canonicalize to lower case, split on nonalphanumerics,
      # and count things that look like matlab identifiers.

      tr/A-Z/a-z/;
      grep($hsymbols{$_}++,grep(/[a-z]\w*/,split('\W',$_)))

    } else {

      # Split off and ignore trailing comments
      # Split on nonalphanumerics and count identifiers

      ($statement,$comment) = split('%',$_,2);
      grep($csymbols{$_}++,grep(/[a-zA-Z]\w*/,split('\W',$statement)));
    }
  }
  close MFILE;

  # Now mark each name that appears in the list of symbols

  # Compute the names appearing among the symbols in the code section
  # (@csym), in the help section (@hsym)

#  if ($opt_g) {   We always want cross listed
  if (1) {

    # If global, search among all symbols

    @csym = grep($csymbols{$_},@names);
    @hsym = grep($hsymbols{$_},@names);

  } else {

    # If not global, search among files appearing in the same directory

    $dir = $mdir{$file};
    @csym = grep(($csymbols{$_} && ($dir eq $mdir{$mfile{$_}})), @names);
    @hsym = grep(($hsymbols{$_} && ($dir eq $mdir{$mfile{$_}})), @names);

  }

  # Now record the cross references in $cref and $href

  grep($cref{$file,$mfile{$_}}=$csymbols{$_},@csym);
  grep($href{$file,$mfile{$_}}=$hsymbols{$_},@hsym);

  undef(%csymbols);
  undef(%hsymbols);
}

######################################################################
#
#  Setup the html directories
#
######################################################################

# Create an html subdirectory name for every unique matlab directory in the
# list @mdirs. The name is constructed using the tail of the directory
# prefaced by a unique number.
#
#  $hdir{$mdir}   - html subdirectory for matlab directory $mdir

$x = 1;
foreach (@mdirs) {
  @z = reverse(split("/",$_));
  $hdir{$_} = "$x.".@z[0];
  $x++;
}

# for each .m file, name a corresponding .html file

foreach (@mfiles) {
  $hfile{$file} = $name{$_}.".html";
}

# Now test a build the corresponding html directories.

print "Checking HTML directories.\n" if $verbose;

if (!-e $hroot) {
  mkdir($hroot,umask) || die("Cannot create directory $hroot\n");
  chmod 0755, $hroot;
}
opendir(HDIR,$hroot) || die ("Cannot open directory $hroot\n");
closedir(HDIR);
die("HTML directory $hroot is not writable\n") if !-w $hroot;

print "HTML Directory $hroot is OK\n" if $verbose;

foreach (@mdirs) {
  local($x) = $hroot."/".$hdir{$_};
  if (!-e $x) {
    mkdir($x,umask) || die("Cannot create directory $x\n");
    chmod(0755,$x);
  }
  opendir(HDIR,$x) || die ("Cannot open directory $x\n");
  closedir(HDIR);
  die("HTML directory $x is not writable\n") if !-w $x;
  print "HTML Directory $x is OK\n" if $verbose;
}

######################################################################
#
#  If verbose mode, write a sparse matrix with the xref data
#
######################################################################

# For kicks, let's write out name and cross reference info in a
# sparse matrix format. These can be anaylzed using matlab.

if ($verbose) {
  local($x,$y,$i,$j);

  open(TXT,">$hroot/names.txt");
  foreach (@mfiles) {
    $b = " " x (20-length($name{$_}));
    print TXT $name{$_},$b,$_,"\n";
  }
  close TXT;

  open(DAT,">$hroot/xref.dat");
  $i = 0;
  foreach $x (@mfiles) {
    $j = 0;
    $i++;
    foreach $y (@mfiles) {
      $j++;
      print DAT "$i  $j  $cref{$x,$y}\n" if $cref{$x,$y};
    }
  }
  close(DAT);
}

######################################################################
#
#  Write the master index file
#
######################################################################

$indexfile = "$hroot/index.html";

print "Writing master $indexfile\n" if $verbose;

open(HFILE,">$indexfile") || die("Cannot open index file $indexfile\n");
print HFILE "<TITLE>Matlab Index</TITLE>\n";
print HFILE '<BODY bgcolor = "FFFFFF">';
print HFILE "\n";
print HFILE "<H1>Matlab Index</H1>\n";
&tagline;

# Print a short introduction

# Print directory listing

print HFILE "<HR><H2>Matlab Directory Indices</H2>\n<pre>\n";
print HFILE "<UL>\n";
foreach $dir (@mdirs) {
  print HFILE "<LI>",&anchor_url("$hdir{$dir}/index.html",$dir),"</LI>\n";
}
print HFILE "</UL>\n";

# Include links to every file that was found

print HFILE "<HR><H2>Identifiers found in these directories</H2>\n";

# We'll do this 4 across in alphabetical order
$i = 1;
foreach (@names) {
# changed  a 15 to 20 to give me more space
  $b = " " x (20 - length($_));
  $html = "$hdir{$mdir{$mfile{$_}}}/$_.html";
  print HFILE &anchor_url($html,$_),$b;
  print HFILE "\n" if (0 == $i%4); # changed a 5 to 4
  $i++;
}

print HFILE "<HR></BODY>\n";
close(HFILE);

######################################################################
#
#  Write an index for each html subdirectory
#
######################################################################

@readme = grep(/readme/i,@mfiles);

foreach $dir (@mdirs) {

  $indexfile = "$hroot/$hdir{$dir}/index.html";

  print "Writing an index file $indexfile\n" if $verbose;

  open(HFILE,">$indexfile") || die("Cannot open index file $indexfile\n");
  print HFILE "<TITLE>Index for Directory $dir</TITLE>\n";
  print HFILE '<BODY bgcolor = "FFFFFF">';
  print HFILE "\n";
  print HFILE &anchor_url("../index.html","[Master Index]"),"\n";
  print HFILE "<H1>Index for $dir</H1>\n";
  &tagline;

  # Now look for a Readme.m file, seemingly a Matlab standard. If there
  # is one, then the help portion is included in the index file.

  foreach $file (@readme) {
    next if !($mdir{$file} eq $dir);
    open(MFILE,$file) || die("Cannot open the file $file");

    # Help Cross Reference information

    undef(@zref);
    foreach (@mfiles) {
      push(@zref,$name{$_}) if $href{$file,$_};
    }

    # Look for the matlab help text block

    $headline = "Readme";
    &writehelpblock;
  }

  # Now write the index catalog for the .m files in this directory

  print HFILE "<HR><H2>Matlab files in this Directory</H2>\n<pre>\n";

  foreach $file (grep($dir eq $mdir{$_},@mfiles)) {
 # changed a 15 to a 20 to give me more space.
    $b = " " x (20 - length($name{$file}));
    $html = $name{$file}.".html";
    print HFILE &anchor_url($html,$name{$file}),"$b$lookfor{$file}\n";
  }
  print HFILE "</pre>\n";
  print HFILE "<HR></BODY>\n";

  close(HFILE);
}

######################################################################
#
#  Write an html file for every m-file
#
######################################################################

# Now write the html file for each matlab file. Need to reread each matlab
# file to find the help text. Note that we can't do this in a single loop
# because we want the back reference information, and also some people
# put the help text before the function declarations.  

# Need a list of mfiles with unique identifiers

@umfiles = values(%mfile);

foreach $file (@mfiles) {

  $h = "$hroot/$hdir{$mdir{$file}}/$name{$file}.html";
  
  print "Writing $h\n" if $verbose;
  # Cross Reference information
  # Find list of names.

  undef(@xref);
  undef(@yref);
  undef(@zref);
  foreach (@umfiles) {
    next if ($name{$_} eq $name{$file});
      push(@xref,$name{$_}) if $cref{$file,$_}; # files we call
      push(@yref,$name{$_}) if $cref{$_,$file}; # files that call us
      push(@zref,$name{$_}) if $href{$file,$_}; # comment lines
    }

  open(MFILE,"<$file") || die("Cannot open $file");
  open(HFILE,">$h") || die("Cannot open $h");

#  print HFILE "<TITLE>$hdir{$file}/$hfile{$file}</TITLE>\n";
  print HFILE "<TITLE>$file</TITLE>\n";
  print HFILE '<BODY bgcolor = "FFFFFF">';
  print HFILE "\n";
  print HFILE &anchor_url("../index.html","[Master Index]"),"\n";
  print HFILE &anchor_url("index.html","[Index for $mdir{$file}]"),"\n";
  print HFILE "<H1>$name{$file}</H1>\n";
  print HFILE "<H2>($mdir{$file}/$name{$file}.m)</H2>\n";

  # If this is a function, then write out the first line as a synposis

  if ($mtype{$file}) {
    print HFILE "<HR><H3>Function Synopsis</H3>\n";
    print HFILE "<pre>$synopsis{$file}</pre>\n";
  }

  # Write the help block

  $headline = "Help text";
  &writehelpblock;

  print HFILE "<HR><H3>Cross-Reference Information</H3>" if (@xref || @yref);
  if (@xref) {
    print HFILE "This $mtype{$file} calls\n";
    print HFILE "<pre><UL>\n";
    foreach (sort @xref) {
      $html = "../$hdir{$mdir{$mfile{$_}}}/$_.html";
      $b = " " x (15 - length($_));
      print HFILE "<LI>",&anchor_url($html,$_),"$b$mfile{$_}</LI>\n";
    }
    print HFILE "</UL></pre>\n";
  }
  if (@yref) {
    print HFILE "This $mtype{$file} is called by\n";
    print HFILE "<pre><UL>\n";
    foreach (sort @yref) {
      $html = "../$hdir{$mdir{$mfile{$_}}}/$_.html";
      $b = " " x (15 - length($_));
      print HFILE "<LI>",&anchor_url($html,$_),"$b$mfile{$_}</LI>\n";
    }
    print HFILE "</UL></pre>\n";
  }

  # Include source text if requested
  # Source text is always requested, replace 1 by $opt_i to revert
  if (1) {
    print HFILE "<HR><H3>Listing of $mtype{$file} $file</H3>\n";
    seek(MFILE,0,0);
    print HFILE "<pre>\n";
    $spaces = " " x 4;
    $red = '<font color = "FF0000">';
    $green = '<font color = "008F00">';
    $endFont = '</font>';
    
    while (<MFILE>) {
      &htmlchar;
      s/\t/$spaces/g;  # convert tabs to 4 whitespace characters.
      if (/^\s*%/) {  # checks to see if we're on a comment line
        foreach $z (@zref) {
          next if /<A.*$z.*A>/;
          $html = "../$hdir{$mdir{$mfile{$z}}}/$z.html";
          s/(\W+)($z)(\W+)/$1.&anchor_url($html,$2).$3/egi;
        }
      } else {
		  foreach $key (@keyWords){ # should select out each keyword
			 s/\b$key\b/$green$key$endFont/g;	# and then highlight it if found
        }
        foreach $x (@xref) {
          next if /<A.*$x.*A>/;
          $html = "../$hdir{$mdir{$mfile{$x}}}/$x.html";
          s/(\W+)($x)(\W+)/$1.&anchor_url($html,$2).$3/eg;
          s/^(\s*)($x)(\W+)/$1.&anchor_url($html,$2).$3/eg;
        }
      }
      s/(%.*)/$red\1$endFont/; # should highlight the comments
      print HFILE $_;
    }
    print HFILE "</pre>\n";
  }
  
  # Print a date stamp

  print HFILE "<HR>\n";
  &tagline;
  print HFILE "</BODY>";
  close(MFILE);
  close(HFILE);
}

######################################################################
#
#  Do some final statistics
#
######################################################################

print $anchor_count," HTML anchors were written.\n" if $verbose;
print $page_count," HTML files were written.\n" if $verbose;
