#!/usr/local/bin/perl
#
# A Perl program to convert "bonehead" matlab to C mex files 
# sutiable for use with MATLAB version 4.X
#
# Call the program as:
#
# m2c m_file_name
#
# Note that the ".m" is left off of the m-file name
#
# Dos users call the program as:
#
# perl m2c m_filename
#
# A new file called m_filename.c is created.
#
# This program converts matlab programs into C "mex" programs.
# It places numerous restrictions on the matlab program to be converted.
# It should only be used if necessary in order to speed up do-loop sections
# in your program. It is not intended as a compiler for an entire matlab
# program.
#
# Don't expect miracles. This program will do about 90% of the grunt work
# of converting matlab to C. However, you need to look at the C code
# for a sanity check. You can run the code through "cc -E" to see how the
# numerous DEFINE statements get expanded.
#
# Additional performance can be attained if all of your array indices
# are integer values. Change the DEFINE statements to remove the
# "(int)" conversions and change all of the index variables from
# double to integer. Also, change the "i=i+1" to "i++" if it is an
# integer and being incremented by one. Finally, don't forget the -O 
# switch when compiling with cmex.
#
# The basic strategy is: 
#
#   1) find the place in your m-file which has the non-vectorizable
#   do-loops
#   2) Cut this part out of the program and put it into a separate
#   matlab function m-file
#   3) Test the program to make sure the function works.
#   4) FORTRAN-ize the function using the rules below. You must
#   expand each vector and matrix operation into x(i,j)-looking
#   commands
#   5) Repeat step 3
#   6) Run M2C on the function file and compile into a mex file
#   7) Repeat step 3 with the mex file
#
# Here are the current rules which need to be followed:
#
#    1) No array operations are allowed. Everything MUST be element-by-element
#    2) Only one statement per line. This includes for, if and end statements.
#    3) No MATLAB function calls allowed. Only the function calls listed in the
#    operators section at the bottom of this script are allowed.
#    4) All arrays MUST be declared explicitly using the ones or zeros command.
#    The only exceptions are the right hand side arguments to the fuction.
#    3) Undeclared scalars are OK. 
#    5) No complex operations are supported.
#    6) Length, ones and zeros operations must be the only operation on the 
#    line. ex: a = length(b) + 1; is illegal. The " + 1" will not be seen.
#    7) Array indices must be simple scalars. X(i+1, 2*j) is NOT ok.
#    
#    Here is a simple matlab function which can be used to make sure
#    M2C is working. It is also an example of what your m-file must look
#    like before attempting conversion with M2C.
#
#    function [ a, b] = testmex(c, d)
#    %
#    % Basic test function for developing my
#    % m2c converter.
#    %
#    
#    if nargin < 2
#       error('Needs more arguments')
#    end
#    
#    [m, n] = size(c);  % Trailing Comment
#    a = zeros(m, n);
#    
#    [mm, nn] = size(d);
#    b = ones( mm,nn);
#    
#    e = ones(mm,1);   % New variable
#    
#    for i=1:2:m
#        for j = i:n
#            a(i,j) = c(i,j) + sin( a(i,j) );
#        end
#    end
#    
#    for i=1:mm
#        for j = i:nn
#            if j < 10
#                b(i,j) = d(i,j) + abs( e(i) );
#            else
#                b(i,j) = d(i,j)^e(i) ;
#            end
#        end
#    end
#
	 
&InitializeStrings;
&InitializeOperators;
print $WelcomeString;
 
$file = shift(@ARGV);
 
open(M_FILE, "< $file.m") || die "Can not open $file.m\n";
open(C_FILE, "> $file.c") || die "Can not open $file.c\n";
select(C_FILE);
$| = 1;
select(STDOUT);
 
# Read in the m-file and store it in a buffer called @source.
# Make some basic formatting changes along the way.

while(<M_FILE>){
 
	chop;
 
	# Skip blank lines
 
	if( $_ =~ /^\s*$/ ){ next;}
 
	$num_lines++;
 
	# Convert comment lines to c-style syntax
 
	if( $_=~/^\s*%/){ $_=~s/^\s*%(.*)$/\/\* $1 \*\//;$comment[$num_lines]=1;}
	else{ $comment[$num_lines] = 0; }
 
	# Strip trailing comments, leading and trailing white space
	$_ =~ s/%.*$//; $_ =~ s/^\s+//; $_ =~ s/\s+$//;
 
	# Strip semicolons
	if( $_ =~ /;\s*\S/ ){ print "$MultipleString \t$_\n\n"; }
	$_ =~ s/;//g;
 
	# Convert control structures
 
	if( ! $comment[$num_lines] ){
 
		$_ =~ s/\belse\b/\}\nelse \{/;
		$_ =~ s/\bfor (.*)$/for\( $1 \)\{/;
		$_ =~ s/\bif\b(.*)$/if\( $1 \)\{/;
		$_ =~ s/\bwhile (.*)$/while\( $1 \)\{/;
		$_ =~ s/\belseif (.*)$/}\nelse if\( $1 \)\{/;
		$_ =~ s/\bend\b/\}/;
		$_ =~ s/\&/\&\&/g;
		$_ =~ s/\|/\|\|/g;
		$_ =~ s/\~/\!/g;
	}
 
	# Convert function names
 
	if( ! $comment[$num_lines] ){
 
		$_ =~ s/\babs\b/fabs/;
		$_ =~ s/\brem\b/fmod/;
		$_ =~ s/\bpi\b/3.14159265358979323846/;
		$_ =~ s/\bnargin\b/nrhs/;
		$_ =~ s/\bnargout\b/nlhs/;
	}
 
	if(s/\berror\s*\(\s*'(.*)'/mexError\(" $1 "/){$TextFunction[$num_lines]=1;}
 
	# Change x^2 to exp(x,2)
 
	if( /\^/ ){ 
 
		s/\s+//g;
		$buf = &FixExponentOperator; 
		$_ = $buf;
	}
 
	$source[$num_lines] = $_;
}
 
if( $source[1] !~ /function/ ){ die $FunctionLineError;  }
 
# Get the input and output variables from function definition
# Function definition can look like :
#        function a = SomeFunction(x) 
#  or    function [a, b, c] = SomeFunction(x, y, z)
 
# Pull out the left hand side variables.

$lhs_substr = $source[1];
$lhs_substr =~ s/function *\[*([^\]=]*).*$/$1/;
$lhs_substr =~ s/,/ /g;
@lhs_vars = split(' ',$lhs_substr);
 
# Pull out the right hand side variables

$rhs_substr = $source[1];
$rhs_substr =~ s/.*\((.*)\)/$1/;
$rhs_substr =~ s/,/ /g;
@rhs_vars = split(' ',$rhs_substr);
 
print "\nLhs Variables = ";
 
foreach $var (@lhs_vars){ print "$var\t"; $lhs_arrays{$var} = 1; }
 
print "\nRhs Variables = ";
 
foreach $var(@rhs_vars){ 
	print "$var\t"; 
	$declared_arrays{$var}=1;
	$rhs_arrays{$var} = 1;
}
 
print "\n";
 
# Print the function definition as a comment in the C program

$source[1] = ("/*").($source[1]).("*/");
$comment[1] = 1;
print C_FILE $source[1], "\n\n";
 
# Search through the code for array declarations with "zeros" or "ones".
 
for( $i=2; $i <= $num_lines; $i++){
 
	if($comment[$i] | $TextFunction[$i]){next;}
 
	if($source[$i] =~ /(ones|zeros)/ ){
 
		$buf = $source[$i];
 
		$buf =~ s/^[\W]+//;
		$buf =~ s/(^[a-zA-Z][a-zA-Z0-9_]*)//;
		$new_array=$1; undef $1;
		$declared_arrays{$new_array} = 1;
 
		if( ! $lhs_arrays{$new_array}){ $temporary_arrays{$new_array} = 1; }
 
		print "Declared Array $new_array\n";
	}
 
}
 
# Search through the matrix definitions to see if there are any 
# capital letters. Exit with a message if there are.
 
@VarsWithCaps = grep( /[A-Z]/, @rhs_vars );
push( @VarsWithCaps, grep( /[A-Z]/, @lhs_vars ) );
push( @VarsWithCaps, grep( /[A-Z]/, keys(%declared_arrays) ) );
 
if( $#VarsWithCaps >= 0 ){
 
	print $CapsWarningString;
 
	foreach $variable (@VarsWithCaps){ print "\t$variable\n"; }
 
	print "\n"; exit;
}

# Search through the code looking for undeclared variables.
 
for( $i=1; $i <= $num_lines; $i++){
 
	if($comment[$i] | $TextFunction[$i]){next;}
 
	$buf = $source[$i];
	$buf =~ s/\W/ /g;                    # Strip Punctuation
	@parts = split(' ', $buf);

	foreach $part ( @parts){

		if( $part =~ /^[a-zA-Z][a-zA-Z0-9_]*$/ ){ 

			if ( ! $declared_arrays{$part} && 
				 ! $scalars{$part}        && 
				 ! $reserved{$part}){

				$scalars{$part} = 1; print "New Scalar Variable $part\n";
			}
		}
	}
}
 
# Put in the standard include files
print C_FILE '#include "cmex.h"' , "\n";
print C_FILE '#include <math.h>' , "\n";
 
# Establish some definitions early on. Array references look like:
# define x(i,j) ( *(mxGetPr(x) + (j-1)*mxGetM(x) + i - 1) )
# Also scan throgh the code replacing the x(i,j) with X(i,j) and 
# the X(i) with X(i,1) and x with X(1,1)
 
print C_FILE $MacroExplanationString;
 
foreach $a (keys %declared_arrays ){
 
	($A=$a) =~ tr/a-z/A-Z/;
	print C_FILE "#define $A(i,j) ( *(${a}Pr+((int)(j)-1)*${a}M+(int)(i)-1) )\n";
 
#   Check each source code line to make sure all matrix references have
#   an index following them. If x is a matrix, y=x is illegal. It needs to
#   look like y=x(1,1).

	for( $i=1; $i <= $num_lines; $i++){
 
		if ( $comment[$i] ){ next; }
 
		if( $source[$i] =~ /\b$a\b\s*[^\(]/ 
			 && $source[$i] !~ /=\s*(ones|zeros|size)/ ){

			print $ScalarString ;
			print "Variable in question: <$a> in line\n$source[$i] \n\n";
		}
 
		$source[$i] =~ s/\b$a\s*\(/$A\(/g;
		$source[$i] =~ s/\b$A\(([^,\)]+)\)/$A\( $1, 1 \)/g;  undef $1;
	}
 
}
 
# Print the entry point

print C_FILE "\nmexFunction(nlhs, plhs, nrhs,prhs)
int nlhs, nrhs;\nMatrix *plhs[], *prhs[];\n{\n\n";
print C_FILE "/* Iteration Variable for Internal Use */\n";
print C_FILE "int iii;\n\n";
 
# declare the scalars

print C_FILE "/* Scalar variables */\n";
print C_FILE "double  "; $dbl_str = "";
foreach $a (keys %scalars ){ $dbl_str .= " $a,"; }
chop($dbl_str);
print C_FILE $dbl_str, ";\n\n";
 
# Declare the arrays

print C_FILE "/* Matrix Variables */\n";
print C_FILE "Matrix "; $mat_str = "";
foreach $a (keys %declared_arrays ){ $mat_str .= " *$a,"; }
chop($mat_str);
print C_FILE $mat_str, ";\n";

print C_FILE "double  "; $pr_str = "";
foreach $a (keys %declared_arrays ){ $pr_str .= " *${a}Pr,"; }
chop($pr_str);
print C_FILE $pr_str, ";\n";
 
print C_FILE "int "; $int_str = "";
foreach $a (keys %declared_arrays ){ $int_str .= " ${a}M,"; }
chop($int_str);
print C_FILE $int_str, ";\n\n";
 
print C_FILE "/* Get Information on Right Hand Side Arguments */\n";
$i = 0; 

foreach $var (@rhs_vars){
 
	print C_FILE "$var = prhs[$i];\n"; $i++;
	print C_FILE $var, "Pr = mxGetPr( $var );\n";
	print C_FILE $var, "M  = mxGetM(  $var );\n\n";
}
 
print C_FILE "\n";
 
# We are finally ready to start translating the actual code segments
 
for( $i=2; $i <= $num_lines; $i++){
 
	if( $comment[$i] ){ print C_FILE $source[$i], "\n"; next; }
 
	$buf = $source[$i];
 
	if( $buf =~ /\bfor\b/ ){   # Parse the "for" loop
 
		print C_FILE "\n";
		$range = $buf;
		$range =~ s/for\( (.*) \){//; $range = $1;
 
		$range =~ s/\s*(\w+)\s*=\s*//; $var = $1;
 
		@range = split(":", $range);
		$min = $range[0]; $max = $range[$#range];
 
		if( $#range == 1){$step = 1;}
 
		else{$step = $range[1];}
 
		$str = "$var = $min; $var <= $max; $var += $step";
		print C_FILE "for( $str ){ \n\n";
 
	}
 
	elsif($buf =~ /\bsize\b/){
 
		undef $1; undef $2;
		$buf =~ s/^\W+(\w+)\W+(\w+)//;
		$m = $1; $n = $2;
		$buf =~ s/\W+size\s*\(\s*(.*)\s*\)//; $a = $1;
		print C_FILE "/* Determine the size of $a */\n";
		print C_FILE "$m = mxGetM( $a );\n";
		print C_FILE "$n = mxGetN( $a );\n\n";
 
	}
 
	elsif($buf =~ /\blength\b/){
 
		$buf =~ s/^\W*(\w+)//; $m = $1;
		$buf =~ s/\W*length\s*\(\s*(.*)\s*\)//; $a = $1;
		print C_FILE "/* Determine the Length of $a*/\n";
		print C_FILE "$m = mxGetM( $a );\n";
		print C_FILE "if( $m <= 1 )\n";
		print C_FILE "    $m = mxGetN( $a );\n\n";
 
	}
 
	elsif($buf =~ /\b(ones|zeros)\b/){
 
		$ones = 1;
 
		if( $buf =~ /zeros/ ){$ones = 0;}
 
		$buf =~ s/^\W*(\w+)//; $a=$1; undef $1; undef $2;
 
		$buf =~ s/\W+(ones|zeros)\W+//;
		$buf =~ s/([^,]*)\,(.*)\)\s*$//;
		$m = $1; $n = $2;
 
		print C_FILE "/* Allocate Memory for $a */\n";
		print C_FILE "$a = mxCreateFull((int)$m,(int)$n,REAL);\n";
		print C_FILE "${a}Pr = mxGetPr( $a );\n";
		print C_FILE "${a}M  = mxGetM(  $a );\n\n";
 
		if( $ones){
 
			print C_FILE "/* Set all Values of $a to 0ne */\n";
			print C_FILE "for(iii = 0; iii < (int)($m*$n); iii++)\n";
			print C_FILE "    *( ${a}Pr + iii) = 1.0; \n\n";
		}
 
	}
 
	elsif($buf =~ /\b(if|elseif|else|while)\b/ ){
	  print C_FILE "\n", $buf, "\n";}
 
	elsif($buf =~ /\}\s*$/ ){print C_FILE $buf, "\n";}
 
	else{print C_FILE $buf, ";\n" ;} 
 
}
 
print C_FILE "\n";
 
# Done with code conversion portion of the program
# Perform necessary clean-up actions

foreach $a (keys %temporary_arrays){ print C_FILE "mxFreeMatrix( $a );\n"; }
 
for( $i=0; $i<=$#lhs_vars; $i++){ print C_FILE "plhs[$i]=$lhs_vars[$i];\n";}
 
print C_FILE "}\n";

# Make everything look nice

print "\n\nIf running MSDOS, you will crash now. Dont worry about it.\n\n";
system("indent $file.c");
system("expand -3 $file.c > m2cTempFile");
system("mv m2cTempFile $file.c");

# This is a list of the operators in matlab that M2C recognizes.
# If the operator is not on this list, either don't use it or
# edit the C program to make it run right.

sub InitializeOperators{
 
	@controls = (    
	"for",
	"while",
	"if",        
	"else",        
	"elseif",        
	"end",        
	"break",        
	"return",        
	"error");
 
	@functions = (    
	"zeros",
	"size",
	"ones",
	"eps",
	"length",
	"pi",
	"nrhs",                # nargin
	"nlhs",                # nargout
	"acos",
	"asin",
	"atan",
	"atan2",
	"ceil",
	"cos",
	"cosh",
	"exp",
	"fabs",                 # abs
	"floor",
	"fmod",                 # rem
	"log",
	"log10",
	"pow",
	"sin",
	"sinh",
	"sqrt",
	"tan",
	"tanh",
	"mexError");            
 
	@operators = (  "*", "-", "/", "+", "^");
	@relational = ( "==", "<", "<=", ">", ">=", "!=", "!");
 
	foreach $control  (@controls) { $reserved{$control }=1; }
	foreach $function (@functions){ $reserved{$function}=1; }
 
}    
 
# Used to change x^2 to pow(x,2). Based on the recursive parser 
# described in the Programming Perl book. Will also change 
# ( (x+3)^2 + (y + z)^3 )^2  to  pow((pow((x+3),2) + pow((y+z),3)),2)  
 
sub FixExponentOperator {
 
	local($left, $right, $op, $output);
	$output="";
 
#   Pull out the first thing on the left which is not a (
 
	$left = &ExtractNextItemExp;
 
#   Get the operator (+, -, =, *, / or ^)
 
	while ($op=&GetOperator){
 
#       Pull out the first thing on the right of the operator which is not a (
 
		$right = &ExtractNextItemExp;
 
		if ($op eq "^"){ $output .= "pow($left,$right)"; $left=""; $right="";}
 
		else{ $output .= " $left $op";  $left=$right; }
	}
 
	$output .= " $right";
	$output;
}
 
sub GetOperator{ s/^([=+*\/\-\^])// && $1; }
 
sub ExtractNextItemExp{
 
	local($val);
 
	if( s/^(\d+\.?\d*|\.\d+)// ){ $val = $1; } # Match a number
 
	elsif( s/^(\w+\([^\)]+\))// ){ $val = $1; } # Match "name(ind)" of array
 
	elsif( s/^(\w+)// ){ $val = $1; }           # Matched "name" of scalar
 
	elsif( s/^(\()// ){                         # Match an open paren. 
 
		$val = &FixExponentOperator;    # Recurse !!
		s/^(\))// || print "$ExpParseError\t$_\n\n";
		$val = "( $val )";        #Restore the stripped ()'s
	}
 
	$val;
}
 
sub InitializeStrings{
 
$WelcomeString = "

		  Welcome to the MATLAB-to-C converter.
 
	 We hope your conversion will be a pleasant one.
 
			   Version 1.2  March 25, 1995
 
   Brought to you by Mike Evanoff aka evanoff@radix.com

   PLEASE send any problems or bug reports to me. Try to
   include the m-file and the resulting c-program.
 
";
 
$MacroExplanationString = "
/* These define statements allow the Matlab arrays to be indexed correctly. Matlab
indexes into its arrays column-wise, not row-wise like C. Also, C is zero-based, 
not one-based like Matlab. These defined macros perform the necessary accessing
of the matlab arrays. */
 
";
 
$CapsWarningString = "

ERROR: The following arrays have capital letters in them !!. This MUST be changed 
in the matlab code, then re-run m2c. Sorry I can't just fix this in the m2c converter.
Maybe next version. The problem is that I define macros for each variable in caps
and if the variable name is in caps, there will be conflicts.
 
";

$FunctionLineError = "

ERROR: The very first line of the matlab m-file MUST contain the function 
definition ie: function [x, y] = do_something( a, b, c).
This INCLUDES comment lines.

";

$ScalarString = "

ERROR: You are refering to an array as a scalar !!
Remember: All NARGIN (input) variables are treated as arrays, not scalars.
Try using x(1) rather than just x for any scalar input variables.

";

$MultipleString = "

ERROR: Sorry, only one matlab statement per line is allowed.
There seems to be more than one statement on this line:
";

$ExpParseError = "

ERROR: There was a problem parsing a command containing the ^ operator.
There seems to be a missing ) or some sort of nesting problem in line:

";

}

