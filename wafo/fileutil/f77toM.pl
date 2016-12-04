#!/usr/local/bin/perl5 -w

# 	f77toM is a Perl 5 script to convert one or more F77 or F90
# files to Matlab M-files.  
#	If you have any questions, please contact Chris Cornuelle at
# bob@math.umn.edu, Minnesota Center for Industrial Mathematics
# at the University of Minnesota.  
#	Dr. C. Cornuelle
#	School of Mathematics/MCIM
#	206 Church Street SE
#	Minneapolis, MN 55455 USA

# WARNING: This program is not up to version 1.0 as of July 1998.
#	There are a number of shortcomings and it is not fully tested.
#	RSVP to bob@math.umn.edu with all bug reports.
#	However, don't tell me about how it will not handle GO TOs in
# your Fortran code.  :^)
#	Note also the location of your perl binary may differ from what
# is on line number 1 here.  Adjust as necessary.
#	All development and testing done on a Pentium II with Linux 2.0.34.
# Your mileage may vary.

# Development Diary (incomplete):
# 01.03.96csc:  This will take one string and replace it with
#		another for all local files ending in .htm or .html.
# 04.28.98csc:  Formerly edit_local_files.  Now attempts a rude
#		conversion of F77 code into M-files.
# 05.07.98csc:  The way to handle file IDs.  Since we cannot just
#		fopen a file with the same fid as its F77 counterpart,
#		an array is required, fid_index.  This must be global. 
# 05.19.98csc: This is version 0.11, and may it be the last.  It
#		will parse all of B Cockburn;s conservation laws code.
# 05.19.98csc: Input from J Treadway, improvements.  V 0.16 now.
# 06.16.98csc: Routine expand_implied_do_loop does just that, up 
#		to 3 indices deep.
# 06.17.98csc: Horrible realization that each variable needs to
#		be tracked in value to handle implied do loops.  Joy.
# 07.09.98csc: Now handles open and close.  V 0.26 now.
# 07.10.98csc: Make fid_index a hash to handle variables.
# 07.13.98csc: Apparent need to find and modify any array or
#		subroutine names that mirror Matlab M-files.  An example
#		is the routine name "input."
# 07.16.98csc: Oh.  Need print, and also formats w/out the format
#		statement.  This means ed_read, ed_print, and ed_write.
# 07.20.98ccs: V 0.33 now.  This pretty much does everything, or
#		at least takes a stab at it.
#####################################################################

# Typical usage stuff here.
if ( 0 > $#ARGV ) {  # Should be one arg, at least.
  usage();
  }  # endif 

###################################################################
# Basic plan is to open files one at a time, parse, edit, generate
# lines of output, print to M-file, close.  Rinse, repeat.

for ( $ifl=0;$ifl<=$#ARGV;$ifl++ ) {  # Loop over arg list files.

  $fname = $ARGV[ $ifl ];
  
# Check it out ...
  stat( $fname );
  if ( !(-f _) || (-z _) || !(-r _) || !(-T _) ) {

    print "$fname has problems for editing, skipping.\n";

    }
  else {

#   Copy Fortran file to backup, open file with original name for input.
    $file_bkp = "";$file_bkp .= $fname."\.bkp";
    system("cp -f $fname $file_bkp");
    print "file: $fname ->\n backup file: $file_bkp ->\n";
    $m_file = $fname;
#   Will this name be a problem?
    $m_name = $m_file;
    $m_name =~ s/\.f$//;
    $is_trouble = mfile_test( $m_name );
#   Rename.
    if ( $is_trouble ) { 
      $m_file = $m_name."ff.m"; 
      print "*** File renamed: was $fname, now $m_file\n";
      }
    elsif ( $m_file =~ m/\.f$/ ) { $m_file =~ s/f$/m/; }
    else { $m_file .= ".m"; }
 
###################################################################
#   Read file in, scanning for labels, managing continuations, etc.
#   This is an initial parsing.
    open( F77FL, "<$file_bkp" );
    $label_num = $do_num = 0;
    $lineno = 0;
    foreach ( <F77FL> ) {

#     The first step has to be to lowercase everything except comments.
      if (!(/^\S+/) && !(/format/)) { 
        $dataline = $_; 
        $_ = lc( $dataline );
        }  # endif

#     Find and modify dangerous mirrors of Matlab M-files.
      if ( !(/^\S+/) && 
           (/(call)\s+(\w+)(\(.+\))?/ ||
            /(function)\s+(\w+)(\(.+\))?/ ||
            /(subroutine)\s+(\w+)(\(.+\))?/) ) {

#       So far look just at routines.
        $dummy = $_;
        $is_trouble = mfile_test( $2 );
        if ($is_trouble) { 
          $text_array[ $lineno ] = 
            "*** Routine renamed: was $2, now ".$2."ff\n";
          print $text_array[ $lineno ];
          $lineno++;
          $_ = "      ".$1." ".$2."ff".$3."\n"; 
          }
        else {
          $_ = $dummy;
          }  # endif

        }  # endif

if ( /^\s+include/ ) { print "INCLUDE: $_"; }
      if ( /^\s{5}\S{1}/ ) {  # Continuations.
        $text_array[ $lineno-1 ] = 
          ed_asterix( $text_array[ $lineno-1 ], $_ );
        }  
#     Look for other interesting tidbits.
      elsif (/^\s+(else)?\s+if\s*\(.+\)\s+(\S+)/) {
#      elsif (/^\s+(else)?\s+if\s*\(.+\)\s+(.+)$/) {

#       An anti-continuation caused by a one-line if clause.  Need to 
#     break the clause out into a new line and tag on an "end" so that
#     it will be in a "standard" form.  Motivation same as for continuation
#     handling here.  See ed_if below.
#    Special case where we are seeing the end of a long single-line if clause.

        if ( $2 ne "then" ) {

          $tmpif = $_;
          $tmpif =~ s/\(\s+/\(/g;
          $tmpif =~ s/\s+\)/\)/g;
          $tmpif =~ s/\.\s+/\./g;
          $tmpif =~ s/\s+\./\./g;
          $tmpif =~ s/\s+/ /g;
          @tmparray = split(' ',$tmpif);
#         The if line.
          $text_array[ $lineno++ ] = "      ".shift(@tmparray)." ".shift(@tmparray)." then\n";
#         The statement.
          $tmpif = join(' ',@tmparray);
          $text_array[ $lineno++ ] = "        ".$tmpif."\n";
#         The endif.
          $text_array[ $lineno++ ] = "      endif\n";

          }
        else {
          $text_array[ $lineno++ ] = $_;
          }  # endif
        }
      elsif ( 0 == $ifl && !(/^\S+/) &&
              ( /\w+\s{1}function\s{1}\w+/ || 
                /\s+subroutine\s{1}\w+/ ) ) {  # First file should be main.
        print "\nERROR - First file in argument list should be the \'main\' program.\n$_\n\n";
        usage();
        exit;
        }
      elsif (/^\s+do (\d+) /) {  # Pull the labels from doloops.
 
#       Note this label needs to match another one somewhere ...
        $do_vals[ $do_num ] = $1;
        $do_list{ $do_vals[ $do_num ] } = $_;
        $do_num++;
        $text_array[ $lineno++ ] = $_;
 
        }
      else {  # We can move on.
        $text_array[ $lineno++ ] = $_;
        }  # endif

      }  # end of foreach
    close( F77FL );
    $lineno--;

###################################################################
#   05.07.98csc: Add a search to map file names-IDs.
    $fid_index{5} = "screen";
    $fid_index{6} = "screen";
    for ( $i=0;$i<=$lineno;$i++ ) {

#     Find and collect the labels.
      if ( $text_array[ $i ] =~ m/^\s{1,4}(\d{1,4})\s+(.+)/ ) {

#       Get the label.
        $label_vals[ $label_num ] = $1;
        $tmplabel = $2;
        $tmplabel =~ s/\s+/ /g;
        @tmparray = split(' ',$tmplabel);

#       IMHO all labels must be to continue, go to, or format.
        if ( !($tmparray[0] =~ m/continue/) &&
             !($tmparray[0] =~ m/go/ && $tmparray[1] =~ m/to/) &&
             !($tmparray[0] =~ m/goto/) &&
             !($tmparray[0] =~ m/format/) ) {

#         Strip off the label.
          $text_array[ $i ] =~ s/^\s+\d+/        /;
#         Then insert a comment to help track gotos.
          $commentary = "%      GOTO ALERT, original F77 nearly: ";
          $commentary .= $label_vals[ $label_num ]." ".$tmplabel."\n";
#         Should be done above ...
          for ( $j=$lineno+1;$j>$i;$j-- ) {
            $text_array[ $j ] = $text_array[ $j-1 ];
            }  # end for on j
          $lineno++;
          $text_array[ $i ] = $commentary;
          $i++;
          }
        else {
#         The hash of labelled lines will be global.
          $label_list{ $label_vals[ $label_num ] } = join(' ',@tmparray);
          }  # endif
        $label_num++;

        }
      elsif ($text_array[$i]=~ m/^\s+((else)?\s+if\s*\(.+?\))\s+(.+?)\s*$/ &&
              $3 ne "then" ) {

#       Extra-long one-line if continuations missed above.
          for ( $j=$lineno+2;$j>$i;$j-- ) {
            $text_array[ $j ] = $text_array[ $j-1 ];
            $text_array[ $j-1 ] = $text_array[ $j-2 ];
            }  # end for on j
          $lineno += 2;
          $text_array[$i] = "    ".$1." then\n";
          $text_array[$i+1] = "      ".$3."\n";
          $text_array[$i+2] = "     endif\n";
          $i += 2;
 
        }  # endif

      }  # end for on i
    $label_num--;

###################################################################
#   Now that we have the text, let's parse.
    open( MFL, ">$m_file" );

    for ( $lines=0;$lines<=$lineno;$lines++ ) {  # This is where the work is done.

      $dataline = "\n";
      $_ = $text_array[ $lines ];

#     Now we have buckets-o-if clauses.
      SWITCH: {
       
#       Handle empty line, of course.
        if (/^\s+$/) { $dataline = $_; last SWITCH; }

#       Comment line trumps all else.
        if (/^\S+/) { ($dataline = $_) =~ s/^\S{1}/\%/; last SWITCH; }

#       Must match read/write to I/O fid and to format label (global).
        if (/^\s+read/) { $dataline = ed_read( $_ ); last SWITCH; }
        if (/^\s+write/) { $dataline = ed_write( $_ ); last SWITCH; }
        if (/^\s+print/) { $dataline = ed_print( $_ ); last SWITCH; }

#       Subprograms.
        if (/^\s+.+\s+function\s+\w+/) { $dataline = ed_function( $_ ); last SWITCH; }
        if (/^\s+subroutine\s+\w+/) { $dataline = ed_subroutine( $_ ); last SWITCH; }
        if (/^\s+call\s+\w+/) { $dataline = ed_call( $_ ); last SWITCH; }

#       Declarations.
        if (/^\s+integer \w+/) { $dataline = ed_integer( $_ ); last SWITCH; }
        if (/^\s+real(\*8)? \w+/) { $dataline = ed_real( $_ ); last SWITCH; }
        if (/^\s+dimension/) { $dataline = ed_dimension( $_ ); last SWITCH; }
        if (/^\s+parameter/) { $dataline = ed_parameter( $_ ); last SWITCH; }
        if (/^\s+data/) { $dataline = ed_data( $_ ); last SWITCH; }
        if (/^\s+common/) { $dataline = ed_common( $_ ); last SWITCH; }

#       Conditionals.
        if (/^\s+if/ || /^\s+else/ || /^\s+endif/) { $dataline = ed_if( $_ ); last SWITCH; }
        if (/^\s+endif/) { $dataline = ed_endif( $_ ); last SWITCH; }

#       Loops.
        if (/^\s+do \d+ /) { $dataline = ed_doloop( $_ ); last SWITCH; }

#       File management.
        if (/^\s+open/) { $dataline = ed_open( $_ ); last SWITCH; }
        if (/^\s+close/) { $dataline = ed_close( $_ ); last SWITCH; }
        if (/^\s+inquire/) { $dataline = ed_inquire( $_ ); last SWITCH; }

#       Miscellaneous constructs.
        if (/^\s+include/) { $dataline = ed_include( $_ ); last SWITCH; }
        if (/^\s+return/) { $dataline = ed_return( $_ ); last SWITCH; }
        if (/^\s+go to \d+/ || 
            /^\s+goto \d+/) { $dataline = ed_goto( $_ ); last SWITCH; }

#       Stuff without much meaning in a Matlab world.
        if (/^\s+implicit/) { $dataline = "\n"; last SWITCH; }
        if (/^\s+end/) { $dataline = "\n"; last SWITCH; }
#       Comment out the stray format lines.
        if (/^\s+\d+\s+format\s*\(/) { $dataline = "% Would be: $_"; last SWITCH; }

#       Parse near end as it occurs in other contexts.
        if (/=/) { $dataline = ed_assignment( $_ ); last SWITCH; }
        if (/^\s+stop\s*\n/) { $dataline = ed_stop( $_ ); last SWITCH; }
        if (/^\s+(\d+)?\s+continue/) { $dataline = ed_continue( $_ ); last SWITCH; }

#       Handle the wacky stuff ungracefully.
#        die "\nERROR: f77toM will not handle this >$_";
#        $dataline = minedit( $_ ); last SWITCH;
        $dataline = commentout( $_ ); last SWITCH;

        }  # end switch block

#     Header stuff.
      if ( 0 == $lines ) {
        $headline = "% Generated by f77toM v0.34 [(c) 1998] from original F77 file: $fname\n\n";
        $dataline = $headline.$dataline;
        }  # endif

#     Output - note that we need to catch the first line regardless.
      if ( $dataline =~ /\S+/ || 0 == $lines ) { print( MFL $dataline ); }      
#      if ( $dataline =~ /\S+/ || 0 == $lines ) { print $dataline; }
 
      }  # end for on lines.

#   Get the last line.
    print( MFL $dataline );
#    print $dataline;

    close( MFL );
    print "  M-file: $m_file\n";

    }  #  endif

  }  #  end for over arg list

print "Done.\n";


######################################################################
# The routine subs ...

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 07.16.98csc:  Initially a clone of ed_write.
sub ed_print {

  my($f77_line) = @_;
  my($outline) = "";
  my($do_flag) = 0;
  my($implied_do) = "";
  my(@arglist) = "";
  my(@fmt_matlab) = "";
  my(@fmt_matvar) = "";
  my(@listidl) = "";
  my(@make_loops) = "";
  my(@mloop) = "";
  my(@array_parenclose,@array_parenopen);
  my(@f77_array,@fmt_array,@fmtpar,@fmtstr);
  my(@idl_inc,@idl_init,@idl_lim,@idl_ndx);
  my(@var_array,@var_type);
  my($Mfid);
  my($argnum,$arr_var,$dummy,$eqflag,$f77fid,$f77fmt,$f77label,$f77vars);
  my($fclose_string,$fmt_edited,$fopen_string,$wr_string,$printfilename);
  my($idl_var,$idum,$ifmt,$ilups,$imat,$inum,$ivar);
  my($jidl,$jmat,$jpar,$lim,$listnum,$n_indices,$nargs,$nchars,$nlines,$npar,$nidl,$nstr);
  my($parencount,$pflag,$space_str);

# Begin.

# Is this a * format or not?  
# Could be a string, could be a string variable, or could be *.
  if ( $f77_line =~ m/^\s+print\s+[\"|\']\((.+?)\)[\"|\']\s*\,(.+)$/ ) {
    $f77fmt = $1;
    $f77vars = $2;
    }
  elsif ( $f77_line =~ m/^\s+print\s+(\w+)\s*\,(.+)$/ ) {
#   Darned if I know what to do here!
    $warning = "ed_print(): WARNING, this print formatting not supported\n  $f77_line Please rewrite with the character variable ".$1." expanded.\n\n";
    print $warning;
    $outline = "% ".$warning.$2."\n";
    return $outline;
    }
  elsif ( $f77_line =~ m/^\s+print\s+\*\s*\,(.+)$/ ) {
    $f77fmt = "";
    $f77vars = $1;
    }
  else {  # Houston, we have a problem.
    die("ed_print(): ERROR - Format wrong in $f77_line");
    }  # endif

# Holdovers.
  $f77fid = "6";
  $Mfid = 1;  # Matlab screen output.

# The rest is variables.
  $f77vars =~ s/\s+//g;

# Extract strings from the format and tag them.
  $fmt_edited = $f77fmt;
  $nstr = 0;
  while ( $fmt_edited =~ m/\'(.*?)\'/ ) {
    $fmtstr[ $nstr ] = $1;
#   Need to escape as needed by Matlab.
    $fmtstr[ $nstr ] =~ s/(\\)/$1$1/g;
    $fmtstr[ $nstr ] =~ s/(%)/$1$1/g;
    $fmt_edited =~ s/(\'.*?\')(\')?/X$nstr$2/; 
    if ( $2 ) { $fmt_edited =~ s/(X$nstr)($2)/$1\,$2/; }
    $nstr++;
    }  # end while

# Extract and label paren'd format bits.
  $npar = 0;
  $fmt_edited =~ s/format\((.+)\)$/$1/;  # Extract formats only.
  while ( $fmt_edited =~ m/\((.+?)\)/ ) {
    $fmtpar[ $npar ] = $1;
    $fmt_edited =~ s/(\(.+?\))/Y$npar/;
    $npar++;
    }  # end while
# Replace label with expanded format.
  for ( $jpar=0;$jpar<$npar;$jpar++ ) {
    $dummy = "";
    $fmt_edited =~ m/(\d+)?Y$jpar/;
    if ( !$1 ) { $lim = 1; }
    else { $lim = $1; }
    for ( $idum=0;$idum<$lim;$idum++ ) { $dummy .= $fmtpar[$jpar].","; }
    chop $dummy;
    $fmt_edited =~ s/(\d+)?Y$jpar/$dummy/;
    }  # end for

# Now tidy up and tag a few remaining things ...
  $fmt_edited =~ s/\s+//g;
  $fmt_edited =~ s/\//N/g;
  $fmt_edited =~ s/N([^\,|^\)])/N\,$1/g;
  $fmt_edited =~ s/N([^\,|^\)])/N\,$1/g;
  $fmt_edited =~ s/([^\,|^\(])N/$1\,N/g;
  $fmt_edited =~ s/([^\,|^\(])N/$1\,N/g;

# Place each output line in an array element and process.
  @fmt_array = split(',',$fmt_edited);
  $imat = 0;
  $fmt_matlab[ $imat ] = $space_str = "";
  $nchars = 0;  # Absolute chars.
  if ( !$f77vars && !f77fmt ) {  # Really ought not to happen.
    die ("ed_print(): ERROR, must have vars with * format here");
    }
  elsif ( $f77vars && !f77fmt ) {
    @var_array = split(',',$f77vars);
    for ( $ivar=0;$ivar<=$#var_array;$ivar++ ) {
#     Assumes the F77 convention, and no chars at all.  Poor policy.
      if ( $var_array[$ivar] =~ m/^[a-ho-z]{1}/ ) { $var_type[$ivar] = "real"; }
      elsif ( $var_array[$ivar] =~ m/^[i-n]{1}/ ) { $var_type[$ivar] = "int"; }
      else { $var_type[$ivar] = "error"; }
      }  # end for on ivar
    }
  else {
    $var_type[0] = "";
    }  # endif

  for ( $ifmt=0;$ifmt<=$#fmt_array;$ifmt++ ) {

    ($space_str,$nchars,$imat,\@fmt_matlab,\@fmt_matvar) = 
      parse_write_fmt($space_str,$nchars,$ifmt,$imat,\@fmt_array,\@fmt_matlab,\@fmt_matvar,\@fmtstr,\@var_type);

    }  # end for on ifmt
# Need a \n at end.
  $fmt_matlab[$imat] .= "\\n";

####################################################################
# Based on the nature of the output variable string, we will make up
# an fprintf for matlab to use, one per line.

# Loop over each line in the output format.
  for ( $jmat=0;$jmat<=$imat;$jmat++ ) {

    $dummy = "fprintf(".$Mfid.",\'".$fmt_matlab[ $jmat ]."\'";
    if ( $fmt_matvar[ $jmat ] ) {
      $dummy .= ",".$fmt_matvar[ $jmat ].");\n";
      }
    else { $dummy .= ");\n"; }  # endif
    $outline .= $dummy;

    }  # end for on jmat

# Characterize the variables, if any.  This involves moving along the 
# print line, parsing as we go, usually by parentheses and commas.
  if ( $f77vars ) {  # We got 'em.

    $nargs = 0;
#   How many?
    @var_array = split(',',$f77vars);
    for ( $ivar=0;$ivar<=$#var_array;$ivar++ ) {

      if ( $var_array[$ivar] =~ m/^\(/ ) {
#       This means we have something special, and we need to paste it
#     back together again.  Probably an implied do loop.
        $idl_var = "";
        $parencount = 0;
        $pflag = 1;
        $eqflag = 0;
#       Run to the end, gathering all between the parens.
        while ( $var_array[$ivar] =~ m/[\(|\)]+/ || $pflag ) {
          @array_parenopen = split('\(',$var_array[$ivar]);
          @array_parenclose = split('\)',$var_array[$ivar]);
          if ( $var_array[$ivar] =~ m/\)$/ ) { $#array_parenclose++; }
          $parencount += $#array_parenclose - $#array_parenopen;
# print "$parencount += $#array_parenclose - $#array_parenopen\n";
          if ( 0 <= $parencount ) { $pflag = 0; }
#         Paste together ...
          if ( $var_array[$ivar] =~ m/\)$/ ) {  # Not array args.
            $idl_var .= $var_array[$ivar]."|";
            }
          elsif ( $eqflag || $var_array[$ivar] =~ m/\=/ ) {
            $idl_var .= $var_array[$ivar]."|";
            $eqflag = 1;
            }
          elsif ( !$eqflag && (-1 > $parencount) ) {  # Non-robust ...
            $idl_var .= $var_array[$ivar].",";
            }
          else {
            $idl_var .= $var_array[$ivar]."|";
# print "HOW DID I GET HERE? -> $var_array[$ivar]\n";
            }
          if ( $ivar > $#var_array ) { die "ERROR: paren mismatch in $f77vars"; }
          $ivar++;
          }  # endwhile
        $idl_var .= $var_array[$ivar];
        $idl_var =~ s/[\||\,]$//;
        $arglist[ $nargs++ ] = $idl_var;

        }
      elsif ( $var_array[$ivar] =~ m/^\w+\(.+?\)$/ ) {  # Array/function.
        $arglist[ $nargs++ ] = $var_array[$ivar];
        }
      elsif ( $var_array[$ivar] =~ m/^\w+\(.+?$/ ) {  # Ditto ...
        $arr_var = "";  
        $parencount = 0;
        $pflag = 1;
#       Run to the end, gathering all between the parens.
        while ( $var_array[$ivar] =~ m/[\(|\)]+/ || $pflag ) {
          @array_parenopen = split('\(',$var_array[$ivar]);
          @array_parenclose = split('\)',$var_array[$ivar]);
          if ( $var_array[$ivar] =~ m/\)$/ ) { $#array_parenclose++; }
          $parencount += $#array_parenclose - $#array_parenopen;
          if ( 0 == $parencount ) { $pflag = 0; }
          $arr_var .= $var_array[$ivar];
          if ( $ivar > $#var_array ) { die "ERROR: paren mismatch in $f77vars"; }
          $ivar++;
          }  # endwhile
        $arr_var .= $var_array[$ivar];
        $arglist[ $nargs++ ] = $arr_var[$ivar];
        }
      elsif ( $var_array[$ivar] =~ m/^\w+$/ ) {  # Simple variable.
        $arglist[ $nargs++ ] = $var_array[$ivar];
        }
      elsif ( $var_array[$ivar] =~ m/[\*|\+|\-|\/]/ ) {  # Statement. 
        $arglist[ $nargs++ ] = $var_array[$ivar];
        }
      else {  # Dull stuff.
        $arglist[ $nargs++ ] = $var_array[$ivar];
        }  # endif

      }  # end for on ivar

#   At this stage the args are space-delimited, array indices are comma-
# delimited, and implied do bits are pipe-delimited.

    }  # endif

# Now compare the arglist to the outline to see where the variables,
# if any, need to go.
  $argnum = 0;
  $listnum = 0;
  while ( $outline =~ m/[\,|X|\'|\(]X[X|(\,\')|\)]/ ) {

    if ( $arglist[ $argnum ] =~ m/\|/ ) {  # Implied do land.

      $nidl = 0;
      $arglist[ $argnum ] =~ s/^\(//;
      $arglist[ $argnum ] =~ s/\)$//;
      @listidl = split('\|',$arglist[ $argnum ]);
      $n_indices = 0;
      for ( $jdl=0;$jdl<=$#listidl;$jdl++ ) {

        if ( $listidl[$jdl] =~ m/(\w+)\=(\S+)/ ) {  # An index variable.

          $idl_ndx[$n_indices] = $1;
          $idl_init[$n_indices] = $2;
          $jdl++;
          $idl_lim[$n_indices] = $listidl[$jdl];
          if ( $listidl[$jdl+1] !=~ m/(\w+)\=(\S+)/ ) {  # Increment.
            $idl_inc[$n_indices] = $listidl[$jdl+1];
            $jdl++;
            }  # endif
#         Construct the necessary Matlab loop(s).
          $mloop[$n_indices] = "     for ".$idl_ndx[$n_indices].
                               " = ".$idl_init[$n_indices];
          if ( $idl_inc[$n_indices] ) { $mloop[$n_indices] .= ":".$idl_inc[$n_indices]; }
          $mloop[$n_indices] .= ":".$idl_lim[$n_indices]."\n";
          $n_indices++;

          }
        else {  # It had better be an output variable name ...
          $outline =~ s/([\'|\(]?)X{1}([(\,\')|\)]?)/$1\,$listidl[ $jdl ]$2/;
          $arglist[ $argnum ] = $listidl[ $jdl ];
          $nidl++;
          $argnum++;
          }  # endif

        }  # end for on jdl

#     Now insert the loop(s) into the string of commands.
      @make_loops = split(';\n',$outline);
      $outline = "";
      for ( $nlines=0;$nlines<=$#make_loops;$nlines++ ) {
        if ( $make_loops[$nlines] =~ m/\($idl_ndx[0]\)/ ) {
#         Note that here we decrement as the IDL reads right-to-left.
          for ( $ilups=$n_indices-1;$ilups>=0;$ilups-- ) {
            if ( 0 < $nlines ) { $make_loops[$nlines-1] .= ";\n"; }
            $make_loops[$nlines] = $mloop[ $ilups ]."       ".$make_loops[$nlines];
            $make_loops[$nlines] .= ";\n     end;\n";
            }  # end for ilups
          }  # endif
        $outline .= $make_loops[$nlines]."\n";
        }  # end for on nlines

      }
    else {  # Reality.
      $outline =~ s/([\'|\(]?\,?)X{1}([(\,\')|\)]?)/$1$arglist[$argnum]$2/;
      $outline =~ s/\,(\,$arglist[$argnum])/$1/;
      $argnum++;
      }  # endif

#   Limit here is arbitrary.
    if ( $argnum > 7 ) { die "ERROR: $argnum > 7"; }

    }  # endwhile

# Open file for appending/creation.
  if ( 1 != $Mfid ) {
    $fopen_string = $Mfid." = fopen(".$printfilename.",\'a\');\n";
    $fclose_string = "fclose(".$Mfid.");\n";
    }
  else {
    $fopen_string = "";
    $fclose_string = "";
    }  # endif
  $outline = join('',$fopen_string,$outline,$fclose_string);

  return $outline;
  }  # end of ed_print


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 07.15.98csc:  Silly.
sub ed_return {


  my($f77_line) = @_;
  my($outline) = "ed_return(): ERROR in f77toM\n";

# Begin.

  $outline = "\n          return;\n";

  return $outline;
  }  # end of ed_return


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 07.15.98csc: Comment out the string and also print to screen.
sub commentout {

  my($f77_line) = @_;
  my($outline) = "commentout(): ERROR in f77toM\n";

# Begin.
  $commentline = "% ".$f77_line;
  $f77_line =~ s/\n//g;
  $f77_line =~ s/\'//g;
  $f77_line =~ s/\"//g;
  $displine = "disp(\'*** What is this (quotes deleted)? ".$f77_line."\');\n";

  $outline = $commentline.$displine;

  return $outline;
  }  # end of commentout


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 07.13.98csc: If string is a Matlab command, return a 1, else
#	a 0.  List is not complete.
sub mfile_test {

  my($testword) = @_;
  my($istrouble) = 0;  # Assume no problems.

# Begin.
  $_ = $testword;

  MSWITCH: {

    if (/^addpath$/) { $istrouble = 1; last MSWITCH; }
    if (/^cd$/) { $istrouble = 1; last MSWITCH; }
    if (/^clear$/) { $istrouble = 1; last MSWITCH; }
    if (/^dbclear$/) { $istrouble = 1; last MSWITCH; }
    if (/^dbcont$/) { $istrouble = 1; last MSWITCH; }
    if (/^dbdown$/) { $istrouble = 1; last MSWITCH; }
    if (/^dbmex$/) { $istrouble = 1; last MSWITCH; }
    if (/^dbquit$/) { $istrouble = 1; last MSWITCH; }
    if (/^dbstack$/) { $istrouble = 1; last MSWITCH; }
    if (/^dbstatus$/) { $istrouble = 1; last MSWITCH; }
    if (/^dbstep$/) { $istrouble = 1; last MSWITCH; }
    if (/^dbstop$/) { $istrouble = 1; last MSWITCH; }
    if (/^dbtype$/) { $istrouble = 1; last MSWITCH; }
    if (/^dbup$/) { $istrouble = 1; last MSWITCH; }
    if (/^debug$/) { $istrouble = 1; last MSWITCH; }
    if (/^demo$/) { $istrouble = 1; last MSWITCH; }
    if (/^diary$/) { $istrouble = 1; last MSWITCH; }
    if (/^disp$/) { $istrouble = 1; last MSWITCH; }
    if (/^edit$/) { $istrouble = 1; last MSWITCH; }
    if (/^getenv$/) { $istrouble = 1; last MSWITCH; }
    if (/^help$/) { $istrouble = 1; last MSWITCH; }
    if (/^inmem$/) { $istrouble = 1; last MSWITCH; }
    if (/^input$/) { $istrouble = 1; last MSWITCH; }
    if (/^load$/) { $istrouble = 1; last MSWITCH; }
    if (/^lookfor$/) { $istrouble = 1; last MSWITCH; }
    if (/^more$/) { $istrouble = 1; last MSWITCH; }
    if (/^path$/) { $istrouble = 1; last MSWITCH; }
    if (/^plot$/) { $istrouble = 1; last MSWITCH; }
    if (/^profile$/) { $istrouble = 1; last MSWITCH; }
    if (/^rmpath$/) { $istrouble = 1; last MSWITCH; }
    if (/^quit$/) { $istrouble = 1; last MSWITCH; }
    if (/^save$/) { $istrouble = 1; last MSWITCH; }
    if (/^title$/) { $istrouble = 1; last MSWITCH; }
    if (/^type$/) { $istrouble = 1; last MSWITCH; }
    if (/^unix$/) { $istrouble = 1; last MSWITCH; }
    if (/^ver$/) { $istrouble = 1; last MSWITCH; }
    if (/^what$/) { $istrouble = 1; last MSWITCH; }
    if (/^which$/) { $istrouble = 1; last MSWITCH; }
    if (/^who$/) { $istrouble = 1; last MSWITCH; }

    }  # end of MSWITCH

  return $istrouble;
  }  # end of mfile_test


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 07.09.98csc: Note of course that the inquire statement is ever
#	used by _real_ Fortran programmers ...
#		In Matlab we will use exist to simulate this.
sub ed_inquire {

  my($f77_line) = @_;
  my($outline) = "ed_inquire(): ERROR in f77toM\n";
  my($file_access) = "sequential";
  my($file_action) = "readwrite";
  my($file_blank) = "null";
  my($file_delim) = "delim";
  my($file_exist) = 0;
  my($file_form) = "";
  my($file_formatted) = "";
  my($file_id) = 0;
  my($file_iolength) = 0;
  my($file_iostat) = 0;
  my($file_name) = "";
  my($file_named) = "";
  my($file_nextrc) = 0;
  my($file_number) = 0;
  my($file_opened) = "";
  my($file_pad) = "";
  my($file_position) = "asis";
  my($file_read) = "";
  my($file_readwrite) = "";
  my($file_recl) = 0;
  my($file_sequential) = "";
  my($file_unformatted) = "";
  my($file_write) = "";
  my(@specifiers);
  my($Mfid,$ispec);

# Begin.
  $f77_line =~ s/\s+//g;

# Let's return the original for safety.
  $orig_line = "% Originally: ".$f77_line;
  $f77_line =~ s/inquire\(//;
  $f77_line =~ s/\)\s*$//;
# These might be dangerous ...
  $f77_line =~ s/\'//g;
  $f77_line =~ s/\"//g;

# Parse for the specifiers.  Note that we do here F77 and F90 ...
  @specifiers = split(',',$f77_line);
  for ( $ispec=0;$ispec<=$#specifiers;$ispec++ ) {

    $_ = $specifiers[$ispec];

    QSWITCH: {

#     File ID.  Note there are two ways to get this.  Ordering important.
      if (/unit=(\w+)/) { $file_id = $1; last QSWITCH; }
      if ( 0 == $ispec && /[^\=]/ && /(\w+)/) { $file_id = $1; last QSWITCH; }
#     File nom.
      if (/file=(\S+)/) { $file_name = $1; last QSWITCH; }

#     Less interesting stuff.
      if (/access=(\S+)/) { $file_access = $1; last QSWITCH; }
      if (/action=(\S+)/) { $file_action = $1; last QSWITCH; }
      if (/blank=(\S+)/) { $file_blank = $1; last QSWITCH; }
      if (/delim=(\S+)/) { $file_delim = $1; last QSWITCH; }
      if (/exist=(\S+)/) { $file_exist = $1; last QSWITCH; }
      if (/form=(\S+)/) { $file_form = $1; last QSWITCH; }
      if (/formatted=(\S+)/) { $file_formatted = $1; last QSWITCH; }
      if (/iolength=(\S+)/) { $file_iolength = $1; last QSWITCH; }
      if (/iostat=(\S+)/) { $file_iostat = $1; last QSWITCH; }
      if (/name=(\S+)/) { $file_name = $1; last QSWITCH; }
      if (/named=(\S+)/) { $file_named = $1; last QSWITCH; }
      if (/nextrc=(\S+)/) { $file_nextrc = $1; last QSWITCH; }
      if (/number=(\S+)/) { $file_number = $1; last QSWITCH; }
      if (/opened=(\S+)/) { $file_opened = $1; last QSWITCH; }
      if (/pad=(\S+)/) { $file_pad = $1; last QSWITCH; }
      if (/position=(\S+)/) { $file_position = $1; last QSWITCH; }
      if (/read=(\S+)/) { $file_read = $1; last QSWITCH; }
      if (/readwrite=(\S+)/) { $file_readwrite = $1; last QSWITCH; }
      if (/recl=(\S+)/) { $file_recl = $1; last QSWITCH; }
      if (/sequential=(\S+)/) { $file_sequential = $1; last QSWITCH; }
      if (/unformatted=(\S+)/) { $file_unformatted = $1; last QSWITCH; }
      if (/write=(\S+)/) { $file_write = $1; last QSWITCH; }

      }  # end of QSWITCH
    
    }  # end for on ispec

# Now create the corresponding Matlab code.
  if ( !$file_name ) { $file_name = $fid_index{ $file_id }; }
  $outline = $orig_line."     inquiry = exist(\'".$file_name."\',\'file\');     exist = inquiry;     name = \'".$file_name."\';\n";

  return $outline;
  }  # end of ed_inquire


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 07.09.98csc: File close.
sub ed_close {

  my($f77_line) = @_;
  my($outline) = "ed_close(): ERROR in f77toM\n";
  my($file_id) = 0;
  my($file_iostat) = 0;
  my($file_status) = "keep";
  my($error_string, $fclose_string) = "";
  my($Mfid,$ispec);

# Begin.
  $f77_line =~ s/\s+//g;
  $f77_line =~ s/close\(//;
  $f77_line =~ s/\)\s*$//;

# Parse for the specifiers.  Note that we do here F77 and F90 ...
  @specifiers = split(',',$f77_line);
  for ( $ispec=0;$ispec<=$#specifiers;$ispec++ ) {

    $_ = $specifiers[$ispec];

    SWITCHRU: {

#     File ID.  Note there are two ways to get this.  Ordering important.
      if (/unit=(\w+)/) { $file_id = $1; last SWITCHRU; }
      if ( 0 == $ispec && /[^\=]/ && /(\w+)/) { $file_id = $1; last SWITCHRU; }
#     Less interesting stuff.
      if (/iostat=(\S+)/) { $file_iostat = $1; last SWITCHRU; }
      if (/status=(\S+)/) { $file_status = $1; last SWITCHRU; }

      }  # end of SWITCHRU

    }  # end for on ispec

# Apply Matlab to infected area.
  $file_name = $fid_index{ $file_id };
  if ( $file_status eq "delete" ) {  # Well, that's different.
    $outline = "     delete(\'".$file_name."\');\n";
    }
  else {
    $Mfid = "fid".$file_id;
    $fclose_string = "     rstatus = fclose(".$Mfid.");\n";
    $error_string = "     if ( 0 > rstatus )\n       error(\'fclose failed with file=".$file_name.", fid=".$Mfid."\');\n     end;\n\n";
    $outline = $fclose_string.$error_string;
    }  # endif

  return $outline;
  }  # end of ed_close


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 07.08.98csc: File open.  Should be straightforward ...
#	open([unit=],iostat=,err=,file=,status=,access=,form=,recl=,blank=)
sub ed_open {

  my($f77_line) = @_;
  my($outline) = "ed_open(): ERROR in f77toM\n";
  my($file_access) = "sequential";
  my($file_action) = "readwrite";
  my($file_blank) = "null";
  my($file_delim) = "delim";
  my($file_err) = 0;
  my($file_form) = "";
  my($file_id) = 0;
  my($file_iostat) = 0;
  my($file_name) = "";
  my($file_pad) = "yes";
  my($file_position) = "asis";
  my($file_recl) = 0;
  my($file_status) = "unknown";
  my($error_string, $fopen_string, $orig_line, $permission) = "";
  my(@specifiers);
  my($Mfid,$ispec);

# Begin.
  $f77_line =~ s/\s+//g;

# Let's return the original for safety.
  $orig_line = "% Originally: ".$f77_line;
  $f77_line =~ s/open\(//;
  $f77_line =~ s/\)\s*$//;
# These might be dangerous ...
  $f77_line =~ s/\'//g;
  $f77_line =~ s/\"//g;

# Parse for the specifiers.  Note that we do here F77 and F90 ...
  @specifiers = split(',',$f77_line);
  for ( $ispec=0;$ispec<=$#specifiers;$ispec++ ) {

    $_ = $specifiers[$ispec];

    IPSWITCH: {

#     File ID.  Note there are two ways to get this.  Ordering important.
      if (/unit=(\w+)/) { $file_id = $1; last IPSWITCH; }
      if ( 0 == $ispec && /[^\=]/ && /(\w+)/) { $file_id = $1; last IPSWITCH; }
#     File nom.
      if (/file=(\S+)/) { $file_name = $1; last IPSWITCH; }

#     Less interesting stuff.
      if (/iostat=(\S+)/) { $file_iostat = $1; last IPSWITCH; }
      if (/status=(\S+)/) { $file_status = $1; last IPSWITCH; }
      if (/access=(\S+)/) { $file_access = $1; last IPSWITCH; }
      if (/form=(\S+)/) { $file_form = $1; last IPSWITCH; }
      if (/recl=(\S+)/) { $file_recl = $1; last IPSWITCH; }
      if (/blank=(\S+)/) { $file_blank = $1; last IPSWITCH; }
      if (/position=(\S+)/) { $file_position = $1; last IPSWITCH; }
      if (/action=(\S+)/) { $file_action = $1; last IPSWITCH; }
      if (/delim=(\S+)/) { $file_delim = $1; last IPSWITCH; }
      if (/pad=(\S+)/) { $file_pad = $1; last IPSWITCH; }
      if (/err=(\S+)/) { $file_err = $1; last IPSWITCH; }

      }  # end of IPSWITCH
    
    }  # end for on ispec

# Consistency check.
  if ( !$file_recl && $file_access eq "direct" ) {
    die "ed_open(): ERROR RECL=$file_recl and ACCESS=$file_access inconsistent";
    }  # endif
  if ( ($file_position eq "rewind" || $file_position eq "append") &&
        $file_access eq "direct" ) {
    die "ed_open(): ERROR POSITION=$file_recl and ACCESS=$file_access inconsistent";
    }  # endif

#######################################################################
# Convert some of the above to generate something that makes sense in
# Matlab.
  $Mfid = "fid".$file_id;
  if ( !$file_form ) {  # Set defaults.
    if ( $file_access eq "direct" ) { $file_form = "unformatted"; }
    else { $file_form = "formatted"; }
    }  # endif
  if ( $file_position eq "rewind" ) {
    $outline = "    frewind(".$Mfid.");\n";
    return $outline;
    }

  if ( $file_status eq "old" ) {
    if ( $file_action eq "readwrite" ) { $permission = "r+"; } 
    elsif ( $file_action eq "read" ) { $permission = "r"; }
    elsif (  $file_action eq "write" ) { $permission = "w"; }
    }
  elsif ( $file_status eq "new" ) {
    if ( $file_status eq "write" ) { $permission = "w"; }
    else { $permission = "w+"; }
    }
  elsif ( $file_status eq "replace" ) {
    $permission = "w+"
    }
  elsif ( $file_status eq "scratch" ) {
#   These will have no name.
    }
  elsif ( $file_status eq "unknown" ) {
    if ( $file_action eq "readwrite" ) { $permission = "r+"; } 
    elsif ( $file_action eq "read" ) { $permission = "r"; }
    elsif (  $file_action eq "write" ) { $permission = "w"; }
    }  # endif

  if ( $file_position eq "append" ) {  # Supersedes the above.
    if ( $file_status eq "old" ) {
      if ( $file_action eq "write" ) { $permission = "a"; }
      else { $permission = "a+"; }
      }
    else { 
      $permission = "a+"; 
      }  # endif
    }  # endif
 
# Store the file ID and name for later.  How to handle integer 
# statements here?  Global here.
  if ( $fid_index{ $file_id } =~ m/\S+/ ) {  # Extant.
    }
  else {  # New.
    $fid_index{ $file_id } = $file_name;
    }  # endif

# Brain dead error checking.
  if ( !$permission ) {  # Whoops.
    die "ed_open(): ERROR Matlab file permission somehow never set";
    }  # endif

# Now create the corresponding Matlab code.
  $fopen_string = "\n     [".$Mfid.",message] = fopen(\'".$file_name."\',\'".$permission."\');\n";
  $error_string = "     if ( 0 > ".$Mfid." )\n       disp(message);\n       error(\'fopen failed with file=".$file_name.", fid=".$Mfid."\');\n     end;\n\n";

  $outline = $orig_line.$fopen_string.$error_string;

  return $outline;
  }  # end of ed_open


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 07.07.98csc: Parses the specially modified F77 format statement
#   and places results into useful arrays for later work.  Used
#   in ed_read.  Call as in
#       ( $space_str, $nchars, $imat, *fmt_matlab, *fmt_matvar ) = parse_read_fmt( $space_str, $nchars, $ifmt, $imat, *fmt_array, *fmtstr, *var_type );
#   $imat is the external loop variable.
sub parse_read_fmt {

  my($space_str,$nchars,$ifmt,$imat,$fmt_array,$fmt_matlab,$fmt_matvar,$fmtstr,$var_type) = @_;
  my($inum,$isp,$lim);

# Begin.

    if ( $fmt_array->[$ifmt] eq "N" ) {  # New line.
      $fmt_matlab->[$imat] .= "\\n";
      $imat++;  # New line, new fscanf wanted.
      $fmt_matlab->[$imat] = "";
      $nchars = 0;
      }
    elsif ( $fmt_array->[$ifmt] =~ /^X(\d+)/ ) {  # Strings.
      $fmt_matlab->[$imat] .= "%s";
      if ( $fmt_matvar->[$imat] ) { $fmt_matvar->[$imat] .= ","; }
      $fmt_matvar->[$imat] .= "\'".$fmtstr->[$1]."\'";
      $nchars += length( $fmtstr->[$1] );
      }
    elsif ( $fmt_array->[$ifmt] =~ /^(\d+)x/ ||
            $fmt_array->[$ifmt] =~ /tr(\d+)/ ) {  # Space insertion.
      $space_str = "";
      for ( $isp=0;$isp<$1;$isp++ ) { $space_str .= " "; }
      $fmt_matlab->[$imat] .= "%s";
      if ( $fmt_matvar->[$imat] ) { $fmt_matvar->[$imat] .= ","; }
      $fmt_matvar->[$imat] .= "\'".$space_str."\'";
      $space_str = "";
      $nchars += $isp-1;
      }
    elsif ( $fmt_array->[$ifmt] =~ /^tl(\d+)/ ) {  # Space deletion.
#     Not supported.
      }
    elsif ( $fmt_array->[$ifmt] =~ /^t(\d+)/ ) {  # Space by absolutes.
      $nspaces = $1 - $nchars;
      $space_str = "";
      for ( $isp=0;$isp<$nspaces;$isp++ ) { $space_str .= " "; }
      $fmt_matlab->[$imat] .= "%s";
      if ( $fmt_matvar->[$imat] ) { $fmt_matvar->[ $imat ] .= ","; }
      $fmt_matvar->[$imat] .= "\'".$space_str."\'";
      $space_str = "";
      $nchars += $nspaces;
      }
    elsif ( $fmt_array->[$ifmt] =~ /^(\d+)?i(\d+)(\.\d+)?/ ) {  # Integers.
      if ( !$1 ) { $lim = 1; }
      else { $lim = $1; }
      for ( $inum=0;$inum<$lim;$inum++ ) {
        $fmt_matlab->[$imat] .= "%$2$3i";
        $fmt_matvar->[$imat] .= "X";
        }  # endif
#      $fmt_matvar->[$imat] .= "X";
      $nchars += $lim * $2;
      }
    elsif ( $fmt_array->[$ifmt] =~ /^(\d+)?([d|e|f|g])(\d+)(\.\d+)?/ ) {
#     Floating point formats.
      if ( !$1 ) { $lim = 1; }
      else { $lim = $1; }
#     Note that 1->1 mapping of formats may NOT be valid.
      for ( $inum=0;$inum<$lim;$inum++ ) {
        $fmt_matlab->[$imat] .= "%$3$4$2";
        $fmt_matvar->[$imat] .= "X";
        }  # endif
      $nchars += $lim * $3;
      }
    elsif ( $fmt_array->[$ifmt] =~ /^a(\d+)/ ) {  # Characters.
      $fmt_matlab->[$imat] .= "%$1s";
      $fmt_matvar->[$imat] .= "X";
      $nchars += $1;
      }
    elsif ( $var_type->[0] ) {  # * format ...
      if ( "real" eq $var_type->[$ifmt] ) {
        }
      }
    else {  # Que?
      print "ERROR: $fmt_array->[$ifmt] unrecognized.\n";
      }  # endif

  return($space_str,$nchars,$imat,$fmt_matlab,$fmt_matvar);

  }  # end of parse_read_fmt


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 07.06.98csc: Parses the specially modified F77 format statement
#	and places results into useful arrays for later work.  Used
#	in ed_write.  Call as in
#		( $space_str, $nchars, $imat, \@fmt_matlab, \@fmt_matvar ) = parse_write_fmt( $space_str, $nchars, $ifmt, $imat, \@fmt_array, \@fmtstr, \@var_type );
#	$imat is the external loop variable.
sub parse_write_fmt {

  my($space_str,$nchars,$ifmt,$imat,$fmt_array,$fmt_matlab,$fmt_matvar,$fmtstr,$var_type) = @_;
  my($inum,$isp,$lim);

# Begin.

    if ( $fmt_array->[$ifmt] eq "N" ) {  # New line.
      $fmt_matlab->[$imat] .= "\\n";
      $imat++;  # New line, new fprintf wanted.
      $fmt_matlab->[$imat] = "";
      $nchars = 0;
      }
    elsif ( $fmt_array->[$ifmt] =~ /^X(\d+)/ ) {  # Strings.
      $fmt_matlab->[$imat] .= "%s";
      if ( $fmt_matvar->[$imat] ) { $fmt_matvar->[$imat] .= ","; }
      $fmt_matvar->[$imat] .= "\'".$fmtstr->[$1]."\'";
      $nchars += length( $fmtstr->[$1] );
      }
    elsif ( $fmt_array->[$ifmt] =~ /^(\d+)x/ ||
            $fmt_array->[$ifmt] =~ /tr(\d+)/ ) {  # Space insertion.
      $space_str = "";
      for ( $isp=0;$isp<$1;$isp++ ) { $space_str .= " "; }
      $fmt_matlab->[$imat] .= "%s";
      if ( $fmt_matvar->[$imat] ) { $fmt_matvar->[$imat] .= ","; }
      $fmt_matvar->[$imat] .= "\'".$space_str."\'";
      $space_str = "";
      $nchars += $isp-1;
      }
    elsif ( $fmt_array->[$ifmt] =~ /^tl(\d+)/ ) {  # Space deletion.
#     Not supported.
      }
    elsif ( $fmt_array->[$ifmt] =~ /^t(\d+)/ ) {  # Space by absolutes.
      $nspaces = $1 - $nchars;
      $space_str = "";
      for ( $isp=0;$isp<$nspaces;$isp++ ) { $space_str .= " "; }
      $fmt_matlab->[$imat] .= "%s";
      if ( $fmt_matvar->[$imat] ) { $fmt_matvar->[ $imat ] .= ","; }
      $fmt_matvar->[$imat] .= "\'".$space_str."\'";
      $space_str = "";
      $nchars += $nspaces;
      }
    elsif ( $fmt_array->[$ifmt] =~ /^(\d+)?i(\d+)(\.\d+)?/ ) {  # Integers.
      if ( !$1 ) { $lim = 1; }
      else { $lim = $1; }
      for ( $inum=0;$inum<$lim;$inum++ ) { 
        $fmt_matlab->[$imat] .= "%$2$3i";
        $fmt_matvar->[$imat] .= "X";
        }  # endif
#      $fmt_matvar->[$imat] .= "X";
      $nchars += $lim * $2;
      } 
    elsif ( $fmt_array->[$ifmt] =~ /^(\d+)?([d|e|f|g])(\d+)(\.\d+)?/ ) {  
#     Floating point formats.
      if ( !$1 ) { $lim = 1; }
      else { $lim = $1; }
#     Note that 1->1 mapping of formats may NOT be valid.
      for ( $inum=0;$inum<$lim;$inum++ ) {
        $fmt_matlab->[$imat] .= "%$3$4$2";
        $fmt_matvar->[$imat] .= "X";
        }  # endif
      $nchars += $lim * $3;
      }
    elsif ( $fmt_array->[$ifmt] =~ /^a(\d+)/ ) {  # Characters.
      $fmt_matlab->[$imat] .= "%$1s";
      $fmt_matvar->[$imat] .= "X";
      $nchars += $1;
      }
    elsif ( $var_type->[0] ) {  # * format ...
      if ( "real" eq $var_type->[$ifmt] ) {
        }
      }
    else {  # Que?
      print "parse_write_fmt(): ERROR; This >$fmt_array->[$ifmt]< unrecognized.\n";
      }  # endif

  return($space_str,$nchars,$imat,$fmt_matlab,$fmt_matvar);

  }  # end of parse_write_fmt


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 06.16.98csc:  Input an implied do loop, output the expanded 
#	form.  Only does up to 3 indices, and has other problems.
sub expand_implied_do_loop {

  my($input_line) = @_;
  my($outline) = "ERROR: in expand_implied_do_loop in f77toM\n;";
  my(@do_init) = (1,1,1);
  my(@do_lim) = (0,0,0);
  my(@do_inc) = (1,1,1);
  my($gstring) = "\n";
  my(@do_var,$i0,$i1,$i2,$init_list,$init_statement,$ndx,$nlst);
  my(@var_args,@var_list,@var_name);

# Begin.

  if ( $input_line =~ m/^\((.+\=+.+)\)/ ) {
 
    $init_statement = $1;
#   From right to left, remove the loops recursively.
    $init_list = $init_statement;
    $ndx = 0;
    while ( $init_list =~ m/(\w+)\=(\w+)\,(\w+)(\,\w+)?$/ ) {
  
      $do_var[$ndx] = $1;
      if ( $3 ) { $do_init[$ndx] = $do_ndx[$ndx] = $2; }
      else { $do_init[$ndx] = $do_ndx[$ndx] = 1; }
      if ( $3 ) { $do_lim[$ndx] = $3; }
      else { $do_lim[$ndx] = 0; }
      if ( $4 ) { $do_inc[$ndx] = $4; }
      else { $do_inc[$ndx] = 1; }
      $ndx++;
      $init_list =~ s/\,\w+\=\w+\,\w+(\,\w+)?$//;
      if ( $init_list =~ m/^\((.+\=+.+)\)/ ) { $init_list = $1; }
  
      }  # end while
  
#   List contents - variables, arrays, etc.
    $nlst = 0;
    while ( $init_list ) {
  
      if ( $init_list =~ m/^(\w+\(.+?\))[\,|^\)|^\w+]?/ ) {  # Array.
        $var_list[$nlst] = $1;
        $var_list_id[$nlst] = "array";
        $var_name[$nlst] = $var_list[$nlst];
        $var_name[$nlst] =~ s/\([\w|\,]+\)//;
        $gstring .= "global ".$var_name[$nlst].";\n";
        $var_args[$nlst] = $var_list[$nlst];
        $var_args[$nlst] =~ s/\w+\((.+)\)/$1/;
        $init_list =~ s/^(\w+\(.+?\))/$2/;
        $init_list =~ s/^\,//;
        }
      elsif ( $init_list =~ m/^(\w+)\,?/ ) {  # Scalar.
#       Is this ever really going to get here?
        $var_list[$nlst] = $1;
        $var_list_id[$nlst] = "scalar";
        $init_list =~ s/^(\w+)(\,)?/$2/;
        }
      else {  # Probably another implied do loop.
        die "do_data(): WARNING, $init_list has unusual list.";
        }  # endif
      $init_list =~ s/^\,//;
      $nlst++;
  
      }  # end while
 
#   Now deconvolve.  Here we will for now assume no more than three
#   array indices are used.
    if ( 3 < $ndx ) { print "ed_data(): WARNING, $ndx indices, more than 3 ...\n"; } 
    $outline = "";
    for ( $i0=$do_init[0];$i0<=$do_lim[0];$i0 += $do_inc[0] ) {
      if ( 1 == $ndx ) {
        $outline .= idl_index_match($nlst,\@do_var,$i0,-1,-1,
         \@var_list_id,\@var_list,\@var_name);
        }
      else {
        for ( $i1=$do_init[1];$i1<=$do_lim[1];$i1 += $do_inc[1] ) {
          if ( 2 == $ndx ) {
            $outline .= idl_index_match($nlst,\@do_var,$i0,$i1,-1,
              \@var_list_id,\@var_list,\@var_name);
            }
          else {
            for ( $i2=$do_init[2];$i2<=$do_lim[2];$i2 += $do_inc[2] ) {
              $outline .= idl_index_match($nlst,\@do_var,$i0,$i1,$i2,
                \@var_list_id,\@var_list,\@var_name);
                }  # end for on i2
              }  # endif
            }  # end for on i1
          }  # endif
        }  # end for on i0
  
#   Tidy it up.
    $outline =~ s/\|$//;
 
    }
  else {  # Not an implied do loop.
    print "expand_implied_do_loop(): WARNING, $input_line is not an implied do loop.\n";
    }  # endif  

  return $outline;
  }  # end of expand_implied_do_loop


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 06.15.98csc: Lots of globals, used to match indices in an
#	implied do loop.  Passed referenced arrays.
sub idl_index_match {

  my(@input_line) = @_;
  my($avar_list_id,$avar_list,$avar_name,$do_var_str);  # Pointers.
  my($do_val0,$do_val1,$do_val2);  
  my($idx,$ilst,$nindices,$nlst,$one_cell);
  my($outline) = "";

# Begin.

  ($nlst,$do_var_str,$do_val0,$do_val1,$do_val2,$avar_list_id,$avar_list,$avar_name) = @input_line;

# Determine number of indices.
if ( 0 > $do_val1 ) { $nindices = 1; }
elsif ( 0 > $do_val2 ) { $nindices = 2; }
else { $nindices = 3; }

# Loop over the implied do loop variables list. 
  for ( $ilst=0;$ilst<$nlst;$ilst++ ) {
    if ( $avar_list_id->[$ilst] eq "array" ) {  # Tricky ...
      $one_cell = $avar_list->[$ilst];
      $one_cell =~ s/^\w+\((.+)\)/$1/;  # Get args.
#     Loop over the indices and replace.
      for ( $idx=0;$idx<$nindices;$idx++ ) {

        if ( $one_cell =~ m/[\(|\,]?($do_var_str->[$idx])[\)|\,]?/ ) {
          if ( 0 == $idx ) { $do_var_val = $do_val0; }
          elsif ( 1 == $idx ) { $do_var_val = $do_val1; }
          else { $do_var_val = $do_val2; }
          $one_cell =~ s/$do_var_str->[$idx]/$do_var_val/;
          if ( 0 == $idx ) { $one_cell = $avar_name->[$ilst]."(".$one_cell.")"; }
          }  # endif

        }  # end for on idx

        $outline .= $one_cell."|";
      }  # endif
    }  # end for on ilst

  return $outline;
  }  # end of idl_index_match 


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 05.29.98csc: Translate F77 builtin functions.
sub f77_functions {

  my($input_line) = @_;
  my($outline) = "f77_functions(): ERROR in f77toM\n";
  my(@botarray,@div_array,@toparray);
  my($arg1,$arg2,$arr1,$arr2,$bot,$dividend,$divisor,$divline,$divsubstr,$dummy);
  my($eql,$fixed,$lhs,$top);
  my($flag) = 1;
  my($flag1,$flag2,$ibot,$idiv,$iop,$ip,$itop);

# Begin.
  $divsubstr = "BOB";
  $zorro = "0";

# A bit of massage.
  $input_line =~ s/\n//g;
  $input_line =~ s/\s+//g;
  if ( $input_line =~ m/^(.+?[^=|^>|^<|^\~])\={1}([^=|^>|^<|^\~].+?)$/ ) {
    $tmp1 = $1; 
    $tmp2 = $2; 
#   Note that spaces flanking = are inserted here ...
    $input_line = $tmp1." = ".$tmp2;
    }  # endif

# Down to business.
  $input_line =~ s/d?sqrt(\(.+\))/sqrt$1/g;
  $input_line =~ s/d?abs(\(.+\))/abs$1/g;
  $input_line =~ s/nint(\(.+\))/round$1/g;
  $input_line =~ s/int(\(.+\))/floor$1/g;
# Compulsiveness.
  $input_line =~ s/(\d+\.)$/$1$zorro/g;
  $input_line =~ s/(\d+\.)([^\d])/$1$zorro$2/g;
  $input_line =~ s/\*\*/\^/g;
# sign is the same.

  if ( $input_line =~ m/\// ) {  # Integer division(s)?

# TEST
    @div_array = split('/',$input_line);
    for ( $idiv=0;$idiv<=$#div_array;$idiv += 2 ) {

      $dividend = $div_array[ $idiv ];
      $dividend =~ s/^\S+\s*\=\s*//;
      $dividend =~ s/^\((.+)\)$/$1/;
      $flag1 = 0;
      $divisor = $div_array[ $idiv + 1 ];
      $divisor =~ s/^\((.+)\)$/$1/;
      $flag2 = 0;
      $_ = $dividend;

      DSWITCH: {
        if (/\.+/ || /^[a-lo-z][\w|\_]*[\(.*?\)]*$/) { last DSWITCH; }
        if (/^[\+|\-]*\d+$/ || /^[i-n][\w|\_]*[\(.*?\)]*$/) { $flag1 = 1; last DSWITCH; }
        if (/[\+|\-|\*]/ || /^[\w|\_]+\(/) { $flag1 = 2; last DSWITCH; }
        if (/\S+/) { $flag1 = 3; last DSWITCH; }
        }  # end of DSWITCH

      $_ = $divisor;
      SWITCHD: {
        if (/\.+/ || /^[a-lo-z][\w|\_]*$/) { last SWITCHD; }
        if (/^[\+|\-]*\d+$/ || /^[i-n][\w|\_]*$/) { $flag2 = 1; last SWITCHD; }
        if (/[\+|\-|\*]/ || /^[\w|\_]+\(/) { $flag2 = 2; last SWITCHD; }
        if (/\S+/) { $flag2 = 3; last SWITCHD; }
        }  # end of DSWITCH

      if ( !$flag1 || !$flag2 ) {  # Skip this.
        }
      elsif ( 1 == $flag1 && 1 == $flag2 ) {  # Integers fo sure.
        $div_array[ $idiv ] = "fix(".$dividend;
        $div_array[ $idiv + 1 ] = $divisor.")";
        $idiv += 2;
        }
      elsif ( 1 == $flag1 ) {
        if ( 2 == $flag2 ) {
          @op_array = split('[\+|\-|\*]',$divisor);
          for ( $iop=0;$iop<=$#op_array;$iop++ ) {  # Test.
            if ( !($op_array[$iop] =~ m/^\d+$/) &&
                 !($op_array[$iop] =~ m/^[i-n][\w|\_]*$/) ) {
              $flag1 = $flag2 = 0;
              $iop = $#op_array+1;
              }  # endif
            }  # end for on iop
          if ( $flag1 && $flag2 ) {
            $div_array[ $idiv ] = "fix(".$dividend;
            $div_array[ $idiv + 1 ] = $divisor.")";
            }  # endif
          }
        else {
print "   f77_functions(): Case not handled; $flag1 $flag2: $dividend|$divisor\n";
          }  # endif
        }
      elsif ( 1 == $flag2 ) {
        if ( 2 == $flag1 ) {
          @op_array = split('[\+|\-|\*]',$dividend);
          for ( $iop=0;$iop<=$#op_array;$iop++ ) {  # Test.
            if ( !($op_array[$iop] =~ m/^\d+$/) &&
                 !($op_array[$iop] =~ m/^[i-n][\w|\_]*$/) ) {
              $flag1 = $flag2 = 0;
              $iop = $#op_array+1;
              }  # endif
            }  # end for on iop
          if ( $flag1 && $flag2 ) {
            $div_array[ $idiv ] = "fix(".$dividend;
            $div_array[ $idiv + 1 ] = $divisor.")";
            }  # endif
          }
        else {
print "   f77_functions(): Case not handled; $flag1 $flag2: $dividend|$divisor\n";
          }  # endif
        }
      elsif ( 2 == $flag1 ) {
        if ( 2 == $flag2 ) {
          @op_array = split('[\+|\-|\*]',$dividend);
          for ( $iop=0;$iop<=$#op_array;$iop++ ) {  # Test.
            if ( !($op_array[$iop] =~ m/^\d+$/) &&
                 !($op_array[$iop] =~ m/^[i-n][\w|\_]*$/) ) {
              $flag1 = $flag2 = 0;
              $iop = $#op_array+1;
              }  # endif
            }  # end for on iop
          @op_array = split('[\+|\-|\*]',$divisor);
          for ( $iop=0;$iop<=$#op_array;$iop++ ) {  # Test.
            if ( !($op_array[$iop] =~ m/^\d+$/) &&
                 !($op_array[$iop] =~ m/^[i-n][\w|\_]*$/) ) {
              $flag1 = $flag2 = 0;
              $iop = $#op_array+1;
              }  # endif
            }  # end for on iop
          if ( $flag1 && $flag2 ) {
            $div_array[ $idiv ] = "fix(".$dividend;
            $div_array[ $idiv + 1 ] = $divisor.")";
            }  # endif
          }
        else {
          $dividend =~ m/^(.*)\((.*)$/;  $arr1 = $1;  $arg1 = $2;
          $divisor =~ /^(.*)\)(.*)$/;  $arr2 = $2;  $arg2 = $1;
          $dummy = "      dummy = ".$arg1."/".$arg2;
#         Recursion alert!
          $dummy = f77_functions( $dummy );
          @dummy_array = split(' ',$dummy);
          @dumber = split('/',$dummy_array[2]);
          @dummy_array = split(' ',$div_array[ $idiv ]);
          $div_array[ $idiv ] = $dummy_array[0]." ".$dummy_array[1]." ".$arr1."(".$dumber[0];
          $div_array[ $idiv + 1 ] = $dumber[1].")".$arr2;
          }  # endif
        }
      elsif ( 3 == $flag1 && 3 == $flag2 ) {
print "   f77_functions(): Case not handled; $flag1 $flag2: $dividend|$divisor\n";
        }
      else {
print "   f77_functions(): Case not possible?; $flag1 $flag2: $dividend|$divisor\n";
        }  # endif      

      }  # end for on idiv
    $input_line = join('/',@div_array);
# print "TEST after: $input_line\n";

    }  # endif

  $outline = $input_line;

  return $outline;
  }  # end of f77_functions


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 05.28.98csc: Handle assignment statements, w/F77 functions.
# 07.16.98csc: Ah, need to handle subleties.  Like integer 
#	division ...
sub ed_assignment {

  my($f77_line) = @_;
  my($outline) = "ed_assignment(): ERROR in f77toM\n";

# Begin.
  $outline = "";
  $outline .= "% Edit of assignment: >".$f77_line;
  $f77_line =~ s/\n//;

# Are there any F77 functions here?  We may need to restate in
# Matlab-ese.
  $f77_line = f77_functions( $f77_line );
  $f77_line .= ";\n";
  $outline .= $f77_line;

  return $outline;
  }  # end of ed_assignment


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 05.19.98csc: From INTEGER to a simple Matlab variable.
sub ed_integer {

  my($f77_line) = @_;
  my(@f77dims);
  my($outline) = "ed_dimension(): ERROR in f77toM\n";

# Begin.

# A redefinition ...
  $f77_line =~ s/\n//g;
  $f77_line =~ s/^\s+integer //;

# Named common blocks are not AFAIK used in Matlab.
  $f77_line =~ s/\/\S+\///;
  $f77_line =~ s/\s+//g;

  $f77_line = do_declare( $f77_line );

# In theory we are now done.
  $outline = $f77_line;

  return "$outline\n";
  }  # end of ed_integer


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 05.19.98csc: From REAL to a simple Matlab variable.
sub ed_real {

  my($f77_line) = @_;
  my(@f77dims);
  my($outline) = "ed_dimension(): ERROR in f77toM\n";

# Begin.

# A redefinition ...
  $f77_line =~ s/\n//g;
  $f77_line =~ s/^\s+real(\*8)? //;

# Named common blocks are not AFAIK used in Matlab.
  $f77_line =~ s/\/\S+\///;
  $f77_line =~ s/\s+//g;

  $f77_line = do_declare( $f77_line );

# In theory we are now done.
  $outline = $f77_line;

  return "$outline\n";
  }  # end of ed_real


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 05.19.98csc: From DIMENSION to just declarations.
sub ed_dimension {

  my($f77_line) = @_;
  my(@f77dims);
  my($outline) = "ed_dimension(): ERROR in f77toM\n";

# Begin.

# A redefinition ...
  $f77_line =~ s/\n//g;
  $f77_line =~ s/^\s+dimension//;

# Named common blocks are not AFAIK used in Matlab.
  $f77_line =~ s/\/\S+\///;
  $f77_line =~ s/\s+//g;
  $f77_line = do_declare( $f77_line );

# In theory we are now done.
  $outline = $f77_line;

  return "$outline\n";
  }  # end of ed_dimension


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 05.19.98csc: Turns a drab Fortran construct into a tasty new
#	Matlab dish!  That is, it turns lists into "declarations."
sub do_declare {

  my($input_line) = @_;
  my(@arrs);
  my($arrays,$iarr,$outline);

# Begin.

# One per line, despite what Matlab says.
  $input_line .= "\;\n";

# Separate out the arrays.
  $input_line = "global ".$input_line;
  $input_line =~ s/\)\,(\w+)/\);\nglobal $1/g;
  $input_line =~ s/(\w+)\,(\w+\()/$1;\nglobal $2/g;
  $arrays = $input_line;
# Destroy non-arrays.
  $arrays =~ s/(\n?global (\w+\,?)+;\n)//g;
  $arrays =~ s/;\n/\|/g;
  $arrays =~ s/global //g;
  $arrays =~ s/;//g;

# Convert any/all arrays.
  @arrs = split('\|',$arrays);
  for ( $iarr=0;$iarr<=$#arrs;$iarr++ ) {
    $arrs[$iarr] = declare_array( $arrs[$iarr] );
    $input_line .= $arrs[$iarr];
    }  # end for on iarr

# Destroy array bounds in globals.
  $input_line =~ s/(global \w+)\(.+\)(;\n)/$1$2/g;

# Each must be on its own line.
  $input_line =~ s/(^|\n)global (\w+)\,(\w+)/$1global $2\;\nglobal $3/g;
  $input_line =~ s/(^|\n)global (\w+)\,(\w+)/$1global $2\;\nglobal $3/g;
  $input_line =~ s/(^|\n)global (\w+)\,(\w+)/$1global $2\;\nglobal $3/g;
  $input_line =~ s/(^|\n)global (\w+)\,(\w+)/$1global $2\;\nglobal $3/g;
 
  $outline = "$input_line\n";
  return $outline;

  }  # end of do_declare


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 05.19.98csc: In Matlab this is just a global assignment, eh?
sub ed_data{

  my($f77_line) = @_;
  my(@f77data,@values,@vars);
  my($outline) = "ed_data(): ERROR in f77toM\n";
  my($data,$idat,$implied_do_flag,$value);
  my($gstring) = "\n";

# Begin.
  $f77_line =~ s/\n//g;

# For general interest ...
  $outline = "\n% Original DATA statement: ".$f77_line."\n";

  $f77_line =~ s/\s+//g;
  $f77_line =~ s/data//;
  $f77_line =~ s/\)$//; 

# Write each out as an assignment/global line pair.
  @f77data = split('\/',$f77_line);

# 06.12.98csc: New approach to nested implied do loops.  Will
#	be a subroutine eventually, to expand on the indices.
  if ( $f77data[0] =~ m/^\((.+\=+.+)\)/ ) {
    $f77data[0] = expand_implied_do_loop( $f77data[0] ); 
    $implied_doflag = 1;
    }
  else {  # Regular stuff.
    $f77data[0] =~ s/\,/\|/g;
    $implied_doflag = 0;
    }  # endif

  @vars = split('\|',$f77data[0]);
  @values = split(',',$f77data[1]);
  if ( $implied_doflag ) { $data = ""; $outline .= $gstring; }
  for ( $idat=0;$idat<=$#vars;$idat++ ) {
    if ( !($implied_doflag) ) { $data = "global ".$vars[$idat].";\n"; }
    $value = $vars[$idat]." = ".$values[$idat].";\n";
    $outline .= $data.$value;
    }  # endfor

  return "$outline\n";
  }  # end of ed_data


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 05.18.98csc: Truly the incantation of the doomed.
sub ed_goto {

  my($f77_line) = @_;
  my($outline) = "ed_goto(): ERROR in f77toM\n";

# Begin.
  $f77_line =~ s/\n//;
  $outline = "disp(\'DANGER: GOTO ALERT >".$f77_line."');\n";

  return $outline;
  }  # end of ed_goto


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 05.18.98csc: Subroutine call becomes an M-file ...
# 05.28.98csc: This probably means we need to do "global" here.
# 07.10.98csc: Um, need to actually execute the subprogram as 
#	well ...
sub ed_call {

  my($f77_line) = @_;
  my($outline) = "ed_call(): ERROR in f77toM\n";
  my($argeth,$subprogram);

# Begin.
  $f77_line =~ s/\n//;
  $f77_line =~ s/\s+//g;
  $f77_line =~ s/^call//;

  $subprogram = $f77_line;
  if ( $f77_line =~ m/(\w+)\((.+)\)/ ) {  # An argument list!
#    $f77_line = "      global ".$1;
    $nom = $1;
    $argeth = $2;
    @argarray = split(',',$argeth);
    for ( $iarg=0;$iarg<=$#argarray;$iarg++ ) {
      if ( !($argarray[$iarg] =~ m/[a-z]+/) ) {
        $argarray[$iarg] = "dummy".$nom.$iarg;
        }  # endif
      }  # end for on iarg
    $argeth = join(',',@argarray);
    }
  else {
#    $f77_line = "      global ".$f77_line;
    $argeth = "dummy_arg";
    }  # endif 
#  $f77_line .= ";\n";
  $f77_line = "";

# Now we actually execute ...
  $f77_line .= "     [".$argeth."] = ".$subprogram.";\n";

  $outline = $f77_line;

  return $outline;
  }  # end of ed_call


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 05.15.98csc: In Matlab this is just a global assignment, eh?
sub ed_parameter {

  my($f77_line) = @_;
  my(@dummy,@f77params,$ipar);
  my($outline) = "ed_parameter(): ERROR in f77toM\n";
  my($param);

# Begin.
  $f77_line =~ s/\s+//g;
  $f77_line =~ s/parameter//;
  $f77_line =~ s/^\(//; 
  $f77_line =~ s/\)$//; 

# Write each out as an assignment/global line pair.
  @f77params = split(',',$f77_line);
  $outline = "";
  for ( $ipar=0;$ipar<=$#f77params;$ipar++ ) {

    $f77params[$ipar] .= ";\n";
    ($param,@dummy) = split('=',$f77params[$ipar]);
    $param = "global ".$param.";\n";
    $outline .= $param.$f77params[$ipar];

    }  # endfor 

  return "$outline\n";
  }  # end of ed_parameter


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 05.11.98csc: Note we assume no gotos ...
# 05.29.98csc: Should match all these labels to do loops, etc.
sub ed_continue {

  my($f77_line) = @_;
  my($f77continue,$f77label);
  my($nlbl);
  my($outline) = "ed_continue(): ERROR in f77toM\n";

# Begin.

# Is there a label?
  $f77_line =~ s/\s+/ /g;
# Note that "<label> continue" is entirely replaced here.
  if ( $f77_line =~ m/\d+/ ) {

    ($f77label,$f77continue) = split(' ',$f77_line);
    $f77continue = $label_list{ $f77label };
#   This ought to be the original scanned label ...
    if ( !($f77continue =~ m/continue/) ) { die "ed_continue(): ERROR - This >$f77continue< should have an F77 continue"; } 
    $outline = "% Original continue: >".$f77_line."\n";

#   Only if this matches a do loop do we append an "end" here.
    for ( $nlbl=0;$nlbl<$do_num;$nlbl++ ) {
      if ( $f77label == $do_vals[$nlbl] ) { 
        $outline .= "% Matches:".$do_list{$do_vals[$nlbl]}."\n     end;\n";
        }  # endif
      }  # end for
        
    }
  else {  # A random continue line, God knows why.
    $outline = "% $f77_line\n"
    }  # endif

  return $outline;
  }  # end of ed_continue


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 05.11.98csc: Turn do into for.
# 07.15.98csc: Note that $array_list{ $arrayname } would be handy
#	for doing comparisons to ensure loop is OK.
sub ed_doloop {

  my($f77_line) = @_;
  my(@f77_array);
  my($f77_orig,$f77continue,$f77label);
  my($warning) = "";
  my($increment);
  my($outline) = "ed_doloop(): ERROR in f77toM\n";

# Begin.
  $f77_orig = $f77_line;
  $f77_line =~ s/\s+/ /g;
  @f77_array = split(' ',$f77_line);
  $f77label = $f77_array[1];
  $f77continue = $label_list{ $f77label };
  if ( !($f77continue =~ m/continue/) ) { die "ed_doloop(): ERROR - This >$f77continue< should be an F77 continue"; }

# Now that we have a do loop, edit to suit.
  splice(@f77_array,0,2,"     for");
  $outline = join(' ',@f77_array);
  $outline =~ s/\,/\:/;

# Now look for (in|de)crement not unity.
  if ( $outline =~ m/\,(\d+)/ ) {
    $increment = $1;
    $outline =~ s/\,\d+//;
    $outline =~ s/:/:$increment:/;
    }  # endif
  $outline =~ m/\w+\=(\w+)\:/;
  if ( !($1 =~ m/\w+/) || $1 < 1 ) {
    $warning = "% ed_doloop(): WARNING, initial loop variable may be too small = $1\n";
    }  # endif
  
# Put a header string on it.
  $f77_orig =~ s/^/\% F77 loop was: \>/;
  $outline = $f77_orig.$outline."\n";

  return $outline;
  }  # end of ed_doloop


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 05.07.98csc: The I in I/O.
# 05.18.98csc: Tackle this in full much later.
# 06.24.98csc: That appears to be now.  Note the many convergences
#	with ed_write, of course ...
# 07.08.98csc: Looks like we will need to rely on external file
#	open and close.
# 07.16.98csc: Modified to handle list-directed input.  Sort of.
sub ed_read {

  my($f77_line) = @_;
  my($outline) = "";
  my($verify_line) = "";
  my(@iflines) = ("","","");
  my($do_flag) = 0;
  my($implied_do) = "";
  my(@arglist) = "";
  my(@fmt_matlab) = "";
  my(@fmt_matvar) = "";
  my(@listidl) = "";
  my(@make_loops) = "";
  my(@mloop) = "";
  my(@array_parenclose,@array_parenopen);
  my(@f77_array,@fmt_array,@fmtpar,@fmtstr);
  my(@idl_inc,@idl_init,@idl_lim,@idl_ndx);
  my(@var_array,@var_type);
  my($Mfid);
  my($argnum,$arr_var,$dummy,$eqflag,$f77fid,$f77fmt,$f77label,$f77vars);
  my($fmt_edited);
  my($idl_var,$idum,$ifmt,$ilups,$imat,$inum,$ivar);
  my($jidl,$jmat,$jpar,$lim,$listnum);
  my($n_indices,$nargs,$nchars,$nlines,$npar,$nidl,$nstr);
  my($parencount,$pflag,$readfilename,$space_str,$wr_string);

# Begin.

# Is this a * format or not?  
# Could be a string, could be a string variable, could be *, could be 
# a label pointing to a format statement.
  if ( $f77_line =~ m/^\s+print\s+[\"|\']\((.+?)\)[\"|\']\s*\,(.+)$/ ) {
    $f77fmt = $1;
    $f77vars = $2;
    $f77fid = "5";
    $Mfid = 1;  # Matlab screen output.
    }
  elsif ( $f77_line =~ m/^\s+print\s+(\w+)\s*\,(.+)$/ ) {
#   Darned if I know what to do here!
    $warning = "ed_read(): WARNING, this read formatting not supported\n  $f77_line Please rewrite with the character variable ".$1." expanded.\n\n";
    print $warning;
    $outline = "% ".$warning.$2."\n";
    return $outline;
    }
  elsif ( $f77_line =~ m/^\s+read\s+\*\s*\,(.+)$/ ) {
    $f77fmt = "";
    $f77vars = $1;
    $f77fid = "5";
    $Mfid = 1;  # Matlab screen output.
    }
  else {  # Must be serious formatting.

    $f77_line =~ s/read\s+\(/read\(/;
    @f77_array = split(' ',$f77_line);

#   Divide and conquer.
    $rd_string = shift( @f77_array );

#   Get the file ID and format label.
    $rd_string =~ s/read\(//;
    $rd_string =~ s/\)//;
    ($f77fid,$f77label) = split(',',$rd_string);

    if ( 5 == $f77fid ) {  # From the screen.
      $Mfid = 1;  # Matlab screen output.
      }
    else {  # To a file.
      $Mfid = "fid".$f77fid;
      $readfilename = $fid_index{ $f77fid };
      }  # endif

#   What if the fid is a variable name?
#   TEST
    if ( $f77fid =~ m/[^\d]/ ) {
      $verify_line = 
        "  clear dummy;\n  if (5==".$f77fid."|6==".$f77fid.")\n    ".$Mfid." = 1;\n  else\n    ".$Mfid." = ".$f77fid.";\n  end;\n";
#     In this case will do input.
      $iflines[0] = "    if (5==".$f77fid."|6==".$f77fid.")\n";
      $iflines[1] = "    else\n";
      $iflines[2] = "    end;\n";
      }  # endif

#   This had bloody well better be a format.  Note it is global.
    if ( $f77label ne "*" ) {
      $f77fmt = $label_list{ $f77label };
      if ( !($f77fmt =~ m/format/) ) { die "ed_read(): ERROR - This >$f77fmt< should have an F77 format"; }
      $f77fmt =~ s/\s+$//;
      }
    else { $f77fmt = ""; } 
  
#   The rest is variables.
    $f77vars = join(' ',@f77_array);
    $f77vars =~ s/\s+//g;

    }  # endif

# Extract strings from the format and tag them.
  $fmt_edited = $f77fmt;
  $nstr = 0;
  while ( $fmt_edited =~ m/\'(.*?)\'/ ) {
    $fmtstr[ $nstr ] = $1;
#   Need to escape as needed by Matlab.
    $fmtstr[ $nstr ] =~ s/(\\)/$1$1/g;
    $fmtstr[ $nstr ] =~ s/(%)/$1$1/g;
    $fmt_edited =~ s/(\'.*?\')(\')?/X$nstr$2/; 
    if ( $2 ) { $fmt_edited =~ s/(X$nstr)($2)/$1\,$2/; }
    $nstr++;
    }  # end while

# Extract and label paren'd format bits.
  $npar = 0;
  $fmt_edited =~ s/format\((.+)\)$/$1/;  # Extract formats only.
  while ( $fmt_edited =~ m/\((.+?)\)/ ) {
    $fmtpar[ $npar ] = $1;
    $fmt_edited =~ s/(\(.+?\))/Y$npar/;
    $npar++;
    }  # end while

# Replace label with expanded format.
  for ( $jpar=0;$jpar<$npar;$jpar++ ) {
    $dummy = "";
    $fmt_edited =~ m/(\d+)?Y$jpar/;
    if ( !$1 ) { $lim = 1; }
    else { $lim = $1; }
    for ( $idum=0;$idum<$lim;$idum++ ) { $dummy .= $fmtpar[$jpar].","; }
    chop $dummy;
    $fmt_edited =~ s/(\d+)?Y$jpar/$dummy/;
    }  # end for

# Now tidy up and tag a few remaining things ...
  $fmt_edited =~ s/\s+//g;
  $fmt_edited =~ s/\//N/g;
  $fmt_edited =~ s/N([^\,|^\)])/N\,$1/g;
  $fmt_edited =~ s/N([^\,|^\)])/N\,$1/g;
  $fmt_edited =~ s/([^\,|^\(])N/$1\,N/g;
  $fmt_edited =~ s/([^\,|^\(])N/$1\,N/g;

# Place each output line in an array element and process.
  @fmt_array = split(',',$fmt_edited);
  $imat = 0;
  $fmt_matlab[ $imat ] = $space_str = "";
  $nchars = 0;  # Absolute chars.
# Really ought not to happen.
  if ( !$f77vars && !f77fmt ) {
    die ("ed_read(): ERROR, must have vars with * format here");
    }
  elsif ( $f77vars ) {
#  elsif ( $f77vars && !f77fmt ) {
    @var_array = split(',',$f77vars);
    for ( $ivar=0;$ivar<=$#var_array;$ivar++ ) {
#     Assumes the F77 convention, and no chars at all.  Poor policy.
      if ( $var_array[$ivar] =~ m/^[a-ho-z]{1}/ ) { $var_type[$ivar] = "real"; }
      elsif ( $var_array[$ivar] =~ m/^[i-n]{1}/ ) { $var_type[$ivar] = "int"; }
      else { $var_type[$ivar] = "error"; }
      }  # end for on ivar
    }
  else {
    $var_type[0] = "";
    }  # endif

# What if we have a * format?  Mock up a real format and feed to 
# the monster below.
  if ( !$f77fmt ) {
    for ( $ivar=0;$ivar<=$#var_array;$ivar++ ) {
      $vsize = length( $var_array[$ivar] );
      $fmt_array[$ivar] = "g".$vsize;
      }  # end for on ivar
    }  # endif

  for ( $ifmt=0;$ifmt<=$#fmt_array;$ifmt++ ) {

    ($space_str,$nchars,$imat,\@fmt_matlab,\@fmt_matvar) =
      parse_read_fmt($space_str,$nchars,$ifmt,$imat,\@fmt_array,\@fmt_matlab,\@fmt_matvar,\@fmtstr,\@var_type);

    }  # end for on ifmt
# Need a \n at end.
  $fmt_matlab[$imat] .= "\\n";

####################################################################
# Based on the nature of the input variable string, we will make up
# an fscanf for matlab to use, one per line.

# Loop over each line in the output format.
  for ( $jmat=0;$jmat<=$imat;$jmat++ ) {

#   When fscanf-ing from fid=1, well, we are just waiting for
# user input from the screen ...
    if ( $Mfid eq "1" ) {
      $dummy = "input(\'>\',\'s\'"
      }
    else {
      $dummy = "fscanf(".$Mfid.",\'".$fmt_matlab[ $jmat ]."\'";
      }  # endif
    if ( $fmt_matvar[ $jmat ] && $Mfid ne "1" ) {
      $dummy .= ",[".($#var_array+1).",1]);\n";
      }
    else { $dummy .= ");\n"; }  # endif
    $outline .= $dummy;

    }  # end for on jmat

# Characterize the variables, if any.  This involves moving along the
# read line, parsing as we go, usually by parentheses and commas.
  if ( $f77vars ) {  # We got 'em.

    $nargs = 0;
#   How many?
    @var_array = split(',',$f77vars);
    for ( $ivar=0;$ivar<=$#var_array;$ivar++ ) {

      if ( $var_array[$ivar] =~ m/^\(/ ) {
#       This means we have something special, and we need to paste it
#     back together again.  Probably an implied do loop.
        $idl_var = "";
        $parencount = 0;
        $pflag = 1;
        $eqflag = 0;
#       Run to the end, gathering all between the parens.
        while ( $var_array[$ivar] =~ m/[\(|\)]+/ || $pflag ) {
          @array_parenopen = split('\(',$var_array[$ivar]);
          @array_parenclose = split('\)',$var_array[$ivar]);
          if ( $var_array[$ivar] =~ m/\)$/ ) { $#array_parenclose++; }
          $parencount += $#array_parenclose - $#array_parenopen;
# print "$parencount += $#array_parenclose - $#array_parenopen\n";
          if ( 0 <= $parencount ) { $pflag = 0; }
#         Paste together ...
          if ( $var_array[$ivar] =~ m/\)$/ ) {  # Not array args.
            $idl_var .= $var_array[$ivar]."|";
            }
          elsif ( $eqflag || $var_array[$ivar] =~ m/\=/ ) {
            $idl_var .= $var_array[$ivar]."|";
            $eqflag = 1;
            }
          elsif ( !$eqflag && (-1 > $parencount) ) {  # Non-robust ...
            $idl_var .= $var_array[$ivar].",";
            }
          else {
            $idl_var .= $var_array[$ivar]."|";
# print "HOW DID I GET HERE? -> $var_array[$ivar]\n";
            }
          if ( $ivar > $#var_array ) { die "ERROR: paren mismatch in $f77vars"; }
          $ivar++;
          }  # endwhile
        $idl_var .= $var_array[$ivar];
        $idl_var =~ s/[\||\,]$//;
        $arglist[ $nargs++ ] = $idl_var;

        }
      elsif ( $var_array[$ivar] =~ m/^\w+\(.+?\)$/ ) {  # Array/function.
        $arglist[ $nargs++ ] = $var_array[$ivar];
        }
      elsif ( $var_array[$ivar] =~ m/^\w+\(.+?$/ ) {  # Ditto ...
        $arr_var = "";
        $parencount = 0;
        $pflag = 1;
#       Run to the end, gathering all between the parens.
        while ( $var_array[$ivar] =~ m/[\(|\)]+/ || $pflag ) {
          @array_parenopen = split('\(',$var_array[$ivar]);
          @array_parenclose = split('\)',$var_array[$ivar]);
          if ( $var_array[$ivar] =~ m/\)$/ ) { $#array_parenclose++; }
          $parencount += $#array_parenclose - $#array_parenopen;
          if ( 0 == $parencount ) { $pflag = 0; }
          $arr_var .= $var_array[$ivar];
          if ( $ivar > $#var_array ) { die "ERROR: paren mismatch in $f77vars"; }
          $ivar++;
          }  # endwhile
        $arr_var .= $var_array[$ivar];
        $arglist[ $nargs++ ] = $arr_var[$ivar];
        }
      elsif ( $var_array[$ivar] =~ m/^\w+$/ ) {  # Simple variable.
        $arglist[ $nargs++ ] = $var_array[$ivar];
        }
      elsif ( $var_array[$ivar] =~ m/[\*|\+|\-|\/]/ ) {  # Statement.
        $arglist[ $nargs++ ] = $var_array[$ivar];
        }
      else {  # Dull stuff.
        $arglist[ $nargs++ ] = $var_array[$ivar];
        }  # endif

      }  # end for on ivar

#   At this stage the args are space-delimited, array indices are comma-
# delimited, and implied do bits are pipe-delimited.

    }  # endif

# Now compare the arglist to the outline to see where the variables,
# if any, need to go.
  $argnum = 0;
  $listnum = 0;
  while ( $outline =~ m/[X|\'|\(]X[X|(\,\')|\)]/ ) {

    if ( $arglist[ $argnum ] =~ m/\|/ ) {  # Implied do land.

      $nidl = 0;
      $arglist[ $argnum ] =~ s/^\(//;
      $arglist[ $argnum ] =~ s/\)$//;
      @listidl = split('\|',$arglist[ $argnum ]);
      $n_indices = 0;
      for ( $jdl=0;$jdl<=$#listidl;$jdl++ ) {

        if ( $listidl[$jdl] =~ m/(\w+)\=(\S+)/ ) {  # An index variable.

          $idl_ndx[$n_indices] = $1;
          $idl_init[$n_indices] = $2;
          $jdl++;
          $idl_lim[$n_indices] = $listidl[$jdl];
          if ( $listidl[$jdl+1] !=~ m/(\w+)\=(\S+)/ ) {  # Increment.
            $idl_inc[$n_indices] = $listidl[$jdl+1];
            $jdl++;
            }  # endif
#         Construct the necessary Matlab loop(s).
          $mloop[$n_indices] = "     for ".$idl_ndx[$n_indices].
                               " = ".$idl_init[$n_indices];
          if ( $idl_inc[$n_indices] ) { $mloop[$n_indices] .= ":".$idl_inc[$n_indices]; }
          $mloop[$n_indices] .= ":".$idl_lim[$n_indices]."\n";
          $n_indices++;

          }
        else {  # It had better be an output variable name ...
          $outline =~ s/([\'|\(]?)X{1}([(\,\')|\)]?)/$1\,$listidl[ $jdl ]$2/;
          $arglist[ $argnum ] = $listidl[ $jdl ];
          $nidl++;
          $argnum++;
          }  # endif

        }  # end for on jdl

#     Now insert the loop(s) into the string of commands.
      @make_loops = split(';\n',$outline);
      $outline = "";
      for ( $nlines=0;$nlines<=$#make_loops;$nlines++ ) {
        if ( $make_loops[$nlines] =~ m/\($idl_ndx[0]\)/ ) {
#         Note that here we decrement as the IDL reads right-to-left.
          for ( $ilups=$n_indices-1;$ilups>=0;$ilups-- ) {
            if ( 0 < $nlines ) { $make_loops[$nlines-1] .= ";\n"; }
            $make_loops[$nlines] = $mloop[ $ilups ]."       ".$make_loops[$nlines]
;
            $make_loops[$nlines] .= ";\n     end;\n";
            }  # end for ilups
          }  # endif
        $outline .= $make_loops[$nlines]."\n";
        }  # end for on nlines

      }
    else {  # Reality.
      $outline =~ s/([\'|\(]?\,?)X{1}([(\,\')|\)]?)/$1$arglist[$argnum]$2/;
      $outline =~ s/\,(\,$arglist[$argnum])/$1/;
#      $outline =~ s/([\'|\(]?)X{1}([(\,\')|\)]?)/$1\,$arglist[ $argnum ]$2/;
      $argnum++;
      }  # endif

#   Limit here is arbitrary.
    if ( $argnum > 7 ) { die "ERROR: $argnum > 7"; }

    }  # endwhile

# Cut and paste.
  if ( !($Mfid =~ /\d+/) ) {  # A variable named unit.
    $input_line = "      dummy = input(\'>\',\'s\');\n";
    $input_line .= "      dumber = dummy;\n";
    $input_line .= "      clear dummy;\n";
    $input_line .= "      dummy = sscanf(dumber,\'%g\');\n";
    $outline = $verify_line.$iflines[0].$input_line.$iflines[1].
               "      dummy = ".$outline.$iflines[2];
    }
  else {
    $input_line = "      dummy = sscanf(dumber,\'%g\');\n";
    $outline = $verify_line."     dumber = ".$outline.$input_line;
    }  # endif

  for ( $ivar=0;$ivar<=$#arglist;$ivar++ ) {
    $outline .= "     ".$arglist[$ivar]." = dummy(".($ivar+1).");\n";
    }  # endfor on ivar
  $outline .= "     clear dummy;\n";

  return $outline;
  }  # end of ed_read


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 05.07.98csc: The fortran include might be mapped to the running
#	of an appropriate M-file of the same name ...
# 07.14.98csc: Vive la France!  A more brute force approach, 
#	where the contents of the include file are dumped to the
#	output.  Ugly, but Matlab5.1 is less finicky this way.
sub ed_include {

  my($f77_line) = @_;
  my($outline) = "";

# Begin.
#  chop $f77_line;
  $f77_line =~ s/include//;
  $f77_line =~ s/\'//g;
# New.
  $f77_line =~ s/\s+//g;
  $fname = $f77_line.".m";
  $outline = "\n% INCLUDE file name: ".$f77_line."; contents follow:\n";
  open( INCL, "<$fname" );
  while (<INCL>) { $outline .= $_; }
  close( INCL );
  $outline .= "% INCLUDE file name: ".$f77_line."; contents end.\n\n";

  return $outline;
  }  # end of ed_include


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 05.06.98csc: Very similar to ed_function, oddly enough.
# 05.28.98csc: Note that it is possible to modify the value(s) of 
#	the argument(s), so these need to be returned.
sub ed_subroutine {

  my($f77_line) = @_;
  my(@f77_array);
  my(@arglist);
  my($global_villagers) = "";
  my($outline) = "ed_subroutine(): ERROR in f77toM\n";
#  local($_);

# Begin.
  $f77_line =~ s/\n//;
  $f77_line =~ s/subroutine//;
  $f77_line =~ s/\s+//g;

  if ( $f77_line =~ m/(\w+)\((.+)\)/ ) {  # An argument list!

    $outline = "\[".$2."\] = ".$f77_line.";\n";
    $global_villagers = $2;
    $global_villagers =~ s/,/;\nglobal /g;
    $global_villagers = "global ".$global_villagers.";\n";

    }
  else {
    $outline = "\[dummy\] = ".$f77_line.";\n";
    }  # endif

  $outline = "      function ".$outline;
# So that something is returned ...
  $outline .= $global_villagers;

  return $outline;
  }  # end of ed_subroutine


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 05.06.98csc: Better than the original, with save, etc.
sub ed_stop {

  my($f77_line) = @_;
  my($outline) = "ed_stop(): ERROR in f77toM\n";

# Begin.

# Note that "stop" is entirely replaced here.
  $outline = "save;\nquit;\n";

  return $outline;
  }  # end of ed_stop


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 05.06.98csc: Common goes to global, we hope ...
sub ed_common {

  my($f77_line) = @_;
  my($outline) = "ed_common(): ERROR in f77toM\n";

# Begin.

# A redefinition ...
  $f77_line =~ s/\n//g;
  $f77_line =~ s/^\s+common //;

# Named common blocks are not AFAIK used in Matlab.
  $f77_line =~ s/\/\S+\///;
  $f77_line =~ s/\s+//g;

  $f77_line = do_declare( $f77_line );

# In theory we are now done.
  $outline = $f77_line;
  return $outline;
  }  # end of ed_common


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 05.15.98csc: Turn f77 array declaration into Matlab.  Up to
#	2-D only.
# 05.21.98csc: Needs a feature to only declare array once.
# 05.27.98csc: Extend to 3-D.
sub declare_array {

  my($f77_line) = @_;
  my($outline) = "";
  my($arrayname);
  my($lbnd1,$lbnd2,$lbnd3,$ubnd1,$ubnd2,$ubnd3);

# Begin.

  $outline .= "\% Original declaration as: $f77_line\n";

# Probably paranoid.
  $f77_line =~ s/\s+//g;

# Rescale any zero-index offsets.
  if ( $f77_line =~ m/(\w+)\((\w+):?(\w+)?,?(\w+)?:?(\w+)?,?(\w+)?:?(\w+)?\)/ ) {

    $arrayname = $1;
    if ( !($array_list{ $arrayname }) ) {  # New array.
      $lbnd1 = $2; $ubnd1 = $3;
      $lbnd2 = $4; $ubnd2 = $5;
      $lbnd3 = $6; $ubnd3 = $7;
      if ( $ubnd1 eq "" ) { $ubnd1 = $lbnd1; $lbnd1 = 1; }
      if ( $lbnd2 ne "" && $ubnd2 eq "" ) { $ubnd2 = $lbnd2; $lbnd2 = 1; }
      if ( $lbnd3 ne "" && $ubnd3 eq "" ) { $ubnd3 = $lbnd3; $lbnd3 = 1; }
      $array_list{ $arrayname } = "$lbnd1..$ubnd1|$lbnd2..$ubnd2|$lbnd3..$ubnd3";
      $outline .= $arrayname." = zeros(".$ubnd1."-".$lbnd1."+1";
      if ( $lbnd3 ne "" ) { 
        $outline .= ",".$ubnd2."-".$lbnd2."+1".",".$ubnd3."-".$lbnd3."+1".");"; 
        }
      elsif ( $lbnd2 ne "" && $lbnd3 eq "" ) { 
        $outline .= ",".$ubnd2."-".$lbnd2."+1".");"; 
        }
      else { 
        $outline .= ",1);"; 
        }  # endif
      
      }
    else {
#      $outline .= "global ".$arrayname.";\n";
      }  # endif

    }
  else {
    $f77_line =~ s/\n//g;
    $outline = "% WARNING: >$f77_line< is not a convertible F77 array!\n";
    }  # endif 

  return $outline;
  }  # end of declare_array


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 04.30.98csc: Handle line continuations.
sub ed_asterix {

  my($firstline,$secondline) = @_;
  my($outline);

# Begin.
  $firstline =~ s/\s*\n//g;
  $secondline =~ s/^\s{5}\S\s+//;

  $outline = join('',$firstline,$secondline);

  return($outline);
  }  # end of ed_asterix


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 04.30.98csc: Convert write to Matlab output - tricky.
# 05.01.98csc: Now that we prescan the labels, they can be mated
#	with the write labels to produce a sane output string.
# 06.22.98csc: Handles most implied do loops.  Yech.
sub ed_write {

  my($f77_line) = @_;
  my($outline) = "";
  my($do_flag) = 0;
  my($implied_do) = "";
  my(@arglist) = "";
  my(@fmt_matlab) = "";
  my(@fmt_matvar) = "";
  my(@listidl) = "";
  my(@make_loops) = "";
  my(@mloop) = "";
  my(@array_parenclose,@array_parenopen);
  my(@f77_array,@fmt_array,@fmtpar,@fmtstr);
  my(@idl_inc,@idl_init,@idl_lim,@idl_ndx);
  my(@var_array,@var_type);
  my($Mfid);
  my($argnum,$arr_var,$dummy,$eqflag,$f77fid,$f77fmt,$f77label,$f77vars);
  my($fclose_string,$fmt_edited,$fopen_string,$wr_string,$writefilename);
  my($idl_var,$idum,$ifmt,$ilups,$imat,$inum,$ivar);
  my($jidl,$jmat,$jpar,$lim,$listnum,$n_indices,$nargs,$nchars,$nlines,$npar,$nidl,$nstr);
  my($parencount,$pflag,$space_str);

# Begin.
  $f77_line =~ s/write\s+\(/write\(/;
  @f77_array = split(' ',$f77_line);

# Divide and conquer.
  $wr_string = shift( @f77_array );
## Not up for unformatted output just yet ...
#  if ( $wr_string =~ m/write\(\S{1,2},\*\)/ ) { die("ed_write(): ERROR - Cannot handle $wr_string"); }

# Get the file ID and format label.
  $wr_string =~ s/write\(//;
  $wr_string =~ s/\)//;
  ($f77fid,$f77label) = split(',',$wr_string);

# This had bloody well better be a format.  Note it is global.
  if ( $f77label ne "*" ) { 
    $f77fmt = $label_list{ $f77label };
    if ( !($f77fmt =~ m/format/) ) { die "ed_write(): ERROR - This >$f77fmt< should have an F77 format"; }
    $f77fmt =~ s/\s+$//;
    }
  else { $f77fmt = ""; }

# The rest is variables.
  $f77vars = join(' ',@f77_array);
  $f77vars =~ s/\s+//g;

  if ( '06' eq $f77fid ) {  # To the screen.
    $Mfid = 1;  # Matlab screen output.
    }
  else {  # To a file.
    $Mfid = "fid".$f77fid;
    $writefilename = $fid_index{ $f77fid };
    }  # endif

# Extract strings from the format and tag them.
  $fmt_edited = $f77fmt;
  $nstr = 0;
  while ( $fmt_edited =~ m/\'(.*?)\'/ ) {
    $fmtstr[ $nstr ] = $1;
#   Need to escape as needed by Matlab.
    $fmtstr[ $nstr ] =~ s/(\\)/$1$1/g;
    $fmtstr[ $nstr ] =~ s/(%)/$1$1/g;
    $fmt_edited =~ s/(\'.*?\')(\')?/X$nstr$2/; 
    if ( $2 ) { $fmt_edited =~ s/(X$nstr)($2)/$1\,$2/; }
    $nstr++;
    }  # end while

# Extract and label paren'd format bits.
  $npar = 0;
  $fmt_edited =~ s/format\((.+)\)$/$1/;  # Extract formats only.
  while ( $fmt_edited =~ m/\((.+?)\)/ ) {
    $fmtpar[ $npar ] = $1;
    $fmt_edited =~ s/(\(.+?\))/Y$npar/;
    $npar++;
    }  # end while
# Replace label with expanded format.
  for ( $jpar=0;$jpar<$npar;$jpar++ ) {
    $dummy = "";
    $fmt_edited =~ m/(\d+)?Y$jpar/;
    if ( !$1 ) { $lim = 1; }
    else { $lim = $1; }
    for ( $idum=0;$idum<$lim;$idum++ ) { $dummy .= $fmtpar[$jpar].","; }
    chop $dummy;
    $fmt_edited =~ s/(\d+)?Y$jpar/$dummy/;
    }  # end for

# Now tidy up and tag a few remaining things ...
  $fmt_edited =~ s/\s+//g;
  $fmt_edited =~ s/\//N/g;
  $fmt_edited =~ s/N([^\,|^\)])/N\,$1/g;
  $fmt_edited =~ s/N([^\,|^\)])/N\,$1/g;
  $fmt_edited =~ s/([^\,|^\(])N/$1\,N/g;
  $fmt_edited =~ s/([^\,|^\(])N/$1\,N/g;

# Place each output line in an array element and process.
  @fmt_array = split(',',$fmt_edited);
  $imat = 0;
  $fmt_matlab[ $imat ] = $space_str = "";
  $nchars = 0;  # Absolute chars.
# Really ought not to happen.
  if ( !$f77vars && !f77fmt ) {
    die ("ed_write(): ERROR, must have vars with * format here");
    }
  elsif ( $f77vars && !f77fmt ) {
    @var_array = split(',',$f77vars);
    for ( $ivar=0;$ivar<=$#var_array;$ivar++ ) {
#     Assumes the F77 convention, and no chars at all.  Poor policy.
      if ( $var_array[$ivar] =~ m/^[a-ho-z]{1}/ ) { $var_type[$ivar] = "real"; }
      elsif ( $var_array[$ivar] =~ m/^[i-n]{1}/ ) { $var_type[$ivar] = "int"; }
      else { $var_type[$ivar] = "error"; }
      }  # end for on ivar
    }
  else {
    $var_type[0] = "";
    }  # endif

  for ( $ifmt=0;$ifmt<=$#fmt_array;$ifmt++ ) {

    ($space_str,$nchars,$imat,\@fmt_matlab,\@fmt_matvar) = 
      parse_write_fmt($space_str,$nchars,$ifmt,$imat,\@fmt_array,\@fmt_matlab,\@fmt_matvar,\@fmtstr,\@var_type);

    }  # end for on ifmt
# Need a \n at end.
  $fmt_matlab[$imat] .= "\\n";

####################################################################
# Based on the nature of the output variable string, we will make up
# an fprintf for matlab to use, one per line.

# Loop over each line in the output format.
  for ( $jmat=0;$jmat<=$imat;$jmat++ ) {

    $dummy = "fprintf(".$Mfid.",\'".$fmt_matlab[ $jmat ]."\'";
    if ( $fmt_matvar[ $jmat ] ) {
      $dummy .= ",".$fmt_matvar[ $jmat ].");\n";
      }
    else { $dummy .= ");\n"; }  # endif
    $outline .= $dummy;

    }  # end for on jmat

# Characterize the variables, if any.  This involves moving along the 
# write line, parsing as we go, usually by parentheses and commas.
  if ( $f77vars ) {  # We got 'em.

    $nargs = 0;
#   How many?
    @var_array = split(',',$f77vars);
    for ( $ivar=0;$ivar<=$#var_array;$ivar++ ) {

      if ( $var_array[$ivar] =~ m/^\(/ ) {
#       This means we have something special, and we need to paste it
#     back together again.  Probably an implied do loop.
        $idl_var = "";
        $parencount = 0;
        $pflag = 1;
        $eqflag = 0;
#       Run to the end, gathering all between the parens.
        while ( $var_array[$ivar] =~ m/[\(|\)]+/ || $pflag ) {
          @array_parenopen = split('\(',$var_array[$ivar]);
          @array_parenclose = split('\)',$var_array[$ivar]);
          if ( $var_array[$ivar] =~ m/\)$/ ) { $#array_parenclose++; }
          $parencount += $#array_parenclose - $#array_parenopen;
# print "$parencount += $#array_parenclose - $#array_parenopen\n";
          if ( 0 <= $parencount ) { $pflag = 0; }
#         Paste together ...
          if ( $var_array[$ivar] =~ m/\)$/ ) {  # Not array args.
            $idl_var .= $var_array[$ivar]."|";
            }
          elsif ( $eqflag || $var_array[$ivar] =~ m/\=/ ) {
            $idl_var .= $var_array[$ivar]."|";
            $eqflag = 1;
            }
          elsif ( !$eqflag && (-1 > $parencount) ) {  # Non-robust ...
            $idl_var .= $var_array[$ivar].",";
            }
          else {
            $idl_var .= $var_array[$ivar]."|";
# print "HOW DID I GET HERE? -> $var_array[$ivar]\n";
            }
          if ( $ivar > $#var_array ) { die "ERROR: paren mismatch in $f77vars"; }
          $ivar++;
          }  # endwhile
        $idl_var .= $var_array[$ivar];
        $idl_var =~ s/[\||\,]$//;
        $arglist[ $nargs++ ] = $idl_var;

        }
      elsif ( $var_array[$ivar] =~ m/^\w+\(.+?\)$/ ) {  # Array/function.
        $arglist[ $nargs++ ] = $var_array[$ivar];
        }
      elsif ( $var_array[$ivar] =~ m/^\w+\(.+?$/ ) {  # Ditto ...
        $arr_var = "";  
        $parencount = 0;
        $pflag = 1;
#       Run to the end, gathering all between the parens.
        while ( $var_array[$ivar] =~ m/[\(|\)]+/ || $pflag ) {
          @array_parenopen = split('\(',$var_array[$ivar]);
          @array_parenclose = split('\)',$var_array[$ivar]);
          if ( $var_array[$ivar] =~ m/\)$/ ) { $#array_parenclose++; }
          $parencount += $#array_parenclose - $#array_parenopen;
          if ( 0 == $parencount ) { $pflag = 0; }
          $arr_var .= $var_array[$ivar];
          if ( $ivar > $#var_array ) { die "ERROR: paren mismatch in $f77vars"; }
          $ivar++;
          }  # endwhile
        $arr_var .= $var_array[$ivar];
        $arglist[ $nargs++ ] = $arr_var[$ivar];
        }
      elsif ( $var_array[$ivar] =~ m/^\w+$/ ) {  # Simple variable.
        $arglist[ $nargs++ ] = $var_array[$ivar];
        }
      elsif ( $var_array[$ivar] =~ m/[\*|\+|\-|\/]/ ) {  # Statement. 
        $arglist[ $nargs++ ] = $var_array[$ivar];
        }
      else {  # Dull stuff.
        $arglist[ $nargs++ ] = $var_array[$ivar];
        }  # endif

      }  # end for on ivar

#   At this stage the args are space-delimited, array indices are comma-
# delimited, and implied do bits are pipe-delimited.

    }  # endif

# Now compare the arglist to the outline to see where the variables,
# if any, need to go.
  $argnum = 0;
  $listnum = 0;
  while ( $outline =~ m/[X|\'|\(]X[X|(\,\')|\)]/ ) {

    if ( $arglist[ $argnum ] =~ m/\|/ ) {  # Implied do land.

      $nidl = 0;
      $arglist[ $argnum ] =~ s/^\(//;
      $arglist[ $argnum ] =~ s/\)$//;
      @listidl = split('\|',$arglist[ $argnum ]);
      $n_indices = 0;
      for ( $jdl=0;$jdl<=$#listidl;$jdl++ ) {

        if ( $listidl[$jdl] =~ m/(\w+)\=(\S+)/ ) {  # An index variable.

          $idl_ndx[$n_indices] = $1;
          $idl_init[$n_indices] = $2;
          $jdl++;
          $idl_lim[$n_indices] = $listidl[$jdl];
          if ( $listidl[$jdl+1] !=~ m/(\w+)\=(\S+)/ ) {  # Increment.
            $idl_inc[$n_indices] = $listidl[$jdl+1];
            $jdl++;
            }  # endif
#         Construct the necessary Matlab loop(s).
          $mloop[$n_indices] = "     for ".$idl_ndx[$n_indices].
                               " = ".$idl_init[$n_indices];
          if ( $idl_inc[$n_indices] ) { $mloop[$n_indices] .= ":".$idl_inc[$n_indices]; }
          $mloop[$n_indices] .= ":".$idl_lim[$n_indices]."\n";
          $n_indices++;

          }
        else {  # It had better be an output variable name ...
          $outline =~ s/([\'|\(]?)X{1}([(\,\')|\)]?)/$1\,$listidl[ $jdl ]$2/;
          $arglist[ $argnum ] = $listidl[ $jdl ];
          $nidl++;
          $argnum++;
          }  # endif

        }  # end for on jdl

#     Now insert the loop(s) into the string of commands.
      @make_loops = split(';\n',$outline);
      $outline = "";
      for ( $nlines=0;$nlines<=$#make_loops;$nlines++ ) {
        if ( $make_loops[$nlines] =~ m/\($idl_ndx[0]\)/ ) {
#         Note that here we decrement as the IDL reads right-to-left.
          for ( $ilups=$n_indices-1;$ilups>=0;$ilups-- ) {
            if ( 0 < $nlines ) { $make_loops[$nlines-1] .= ";\n"; }
            $make_loops[$nlines] = $mloop[ $ilups ]."       ".$make_loops[$nlines];
            $make_loops[$nlines] .= ";\n     end;\n";
            }  # end for ilups
          }  # endif
        $outline .= $make_loops[$nlines]."\n";
        }  # end for on nlines

      }
    else {  # Reality.
      $outline =~ s/([\'|\(]?\,?)X{1}([(\,\')|\)]?)/$1$arglist[$argnum]$2/;
      $outline =~ s/\,(\,$arglist[$argnum])/$1/;
#      $outline =~ s/([\'|\(]?)X{1}([(\,\')|\)]?)/$1\,$arglist[ $argnum ]$2/;
      $argnum++;
      }  # endif

#   Limit here is arbitrary.
    if ( $argnum > 7 ) { die "ERROR: $argnum > 7"; }

    }  # endwhile

# Open file for appending/creation.
  if ( 1 != $Mfid ) {
    $fopen_string = $Mfid." = fopen(".$writefilename.",\'a\');\n";
    $fclose_string = "fclose(".$Mfid.");\n";
    }
  else {
    $fopen_string = "";
    $fclose_string = "";
    }  # endif
  $outline = join('',$fopen_string,$outline,$fclose_string);

  return $outline;
  }  # end of ed_write

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 04.30.98csc: endif -> end;
sub ed_endif {

  my($f77_line) = @_;
  my($outline) = "";

# Begin.
  ($outline = $f77_line) =~ s/endif/end;/; 

  return $outline;
  }  # end of ed_endif


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 04.29.98csc: Go to Matlab if ...
# 04.30.98csc: Leading spaces are important for line continuation.
# 07.09.98csc: Add Arithmetic If ...
sub ed_if {

  my($f77_line) = @_;
  my(@f77_array);
  my($outline) = "";
  my($condit,$i);

# Begin.

# Destroy all space.
  $f77_line =~ s/\(\s+/\(/g;
  $f77_line =~ s/\s+\)/\)/g;
  $f77_line =~ s/\.\s+/\./g;
  $f77_line =~ s/\s+\./\./g;
  $f77_line =~ s/(if)(\(.+)/$1 $2/;
  $f77_line =~ s/(.+\))(.+)/$1 $2/;
  @f77_array = split(' ',$f77_line);

  $i=0;
  foreach $_ ( @f77_array ) {

    if ( /elseif/ ) {  # Then the next element may be conditional.

      $outline .= "     elseif ";
      $condit = $f77_array[$i+1];
#     Yeah, pretty ugly.
      $condit =~ s/\.eq\./\=\=/g;
      $condit =~ s/\.ne\./\~\=/g;
      $condit =~ s/\.ge\./\>\=/g;
      $condit =~ s/\.le\./\<\=/g;
      $condit =~ s/\.gt\./\>/g;
      $condit =~ s/\.lt\./\</g;
      $condit =~ s/\.or\./\|/g;
      $condit =~ s/\.and\./\&/g;
      $condit = f77_functions( $condit );
      $outline .= $condit;

      return "$outline\n";

      }
    elsif ( /endif/ ) {
      $outline .= "% endif here\n     end;";
      return "$outline\n";
      } 
    elsif ( /if/ ) {  # Then the next element is the conditional.

      $outline = "% if clause begins here.\n";
      $outline .= "     if ";
      $condit = $f77_array[$i+1];
#     Is this an arithmetic if?
      if ( $condit =~ m/(\d+)\,(\d+)\,(\d+)/ ) {
        $outline = "disp(\'DANGER: ARITHMETIC IF ALERT >".$f77_line."');\n";
        return $outline;
        }  # endif
#     Yeah, pretty ugly.
      $condit =~ s/\.eq\./\=\=/g;
      $condit =~ s/\.ne\./\~\=/g;
      $condit =~ s/\.ge\./\>\=/g;
      $condit =~ s/\.le\./\<\=/g;
      $condit =~ s/\.gt\./\>/g;
      $condit =~ s/\.lt\./\</g;
      $condit =~ s/\.or\./\|/g;
      $condit =~ s/\.and\./\&/g;
      $condit = f77_functions( $condit );
      $outline .= $condit;
      return "$outline\n";

      }
    elsif ( /else/ ) {
      $outline .= "     else ";
      }
    elsif ( /then/ ) {
      }
    else {  # Uh oh.
      die("ed_if(): Why me?  $_");
      }  # endif

    $i++;
    }  # end foreach

# Should not get here, ja?
  return "$outline\n";
  }  # end of ed_if

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 04.29.98csc: Go to Matlab function.
sub ed_function {

  my($f77_line) = @_;
  my(@f77_array);
  my($global_villagers) = "";
  my($outline) = "ERROR: in ed_function in f77toM\n";
  local($_);

# Begin.
  @f77_array = split(' ',$f77_line);
  foreach $_ ( @f77_array ) {

    if ( /real/ || /function/ ) { }
    else {  # Must be the function name and args.

      $_ =~ s/^/function \[dummy\] \= /;
      ($outline = $_) =~ s/$/\;\n/;
      if ( $outline =~ m/\((.+)\)/ ) {

        $global_villagers = $1;
        $outline =~ s/dummy/$global_villagers/;
#        $outline =~ s/dummy/dummy,$global_villagers/;
        $global_villagers =~ s/,/;\nglobal /g;
        $global_villagers = "global ".$global_villagers.";\n";

        }  # endif

      }  # endif

    }  # end foreach

# So that something is returned ...
  $outline .= "global dummy;\n";
  $outline .= $global_villagers;

  return( $outline );
  }  # end of ed_function


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 04.29.98csc: To F77 just add a ; at end and return.
# 05.28.98csc: Do a comment line as well ...
sub minedit {

  my($f77_line) = @_;
  my($outline) = "ERROR: in minedit in f77toM\n";

# Begin.
  $outline = "";
  $outline .= "% Min edit of: >".$f77_line;
  $f77_line =~ s/\s+//g;
  $f77_line .= ";\n";
  $outline .= $f77_line;
  
  return $outline;
  }  # end of minedit


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# 04.29.98csc: Give a usage.
sub usage {

# Begin.
  print "\nUsage: f77toM <list of Fortran files>\n\n";
  print "Takes input Fortran program files, makes a crude conversion \n";
  print "to Matlab, outputs the corresponding M-files.  For every input\n";
  print "file there is an output file with the .m suffix, and a backup.\n";
  print "\nExample:\n  f77toM dog.f cat.f  produces  dog.m cat.m, and dog.f.bkp cat.f.bkp\n"; 
  print "\nIn the argument list, all f77 files that open external files must \n";
  print "occur before f77 files that do I/O to those files.  Generally, use\n";
  print "the ordering: main, includes, openers, others.\n";
  print "\n07.20.98csc, v0.34\n\n";
  exit(2);

  }  # end of usage

1;
