# Makefile.PL for a package defined by PP code.

use PDL::Core::Dev;            # Pick up development utilities
use ExtUtils::MakeMaker;
use Config;
use ExtUtils::F77;

# -- Add new subroutines here! --

my @src = qw(grdlist.c grd2pdl.c);

# -- Try not to mess with below stuff --

my $compiler_available = ExtUtils::F77->testcompiler;
my $myf77    = ExtUtils::F77->compiler();

my $underscore = '';
my $trail_ = ExtUtils::F77->trail_;  # fortran names need trailing underscore
if (($myf77 eq 'g77') and $trail_) { $underscore = '-fno-underscoring'; }
if (($myf77 eq 'f77') and $trail_ and ($Config{'osname'} eq 'dec_osf')) { $underscore = '-assume nounderscore'; }
print "Underscoring define = $underscore\n";
# Might want to add 'no underscore' rules for other architectures below...

my $runtime = ExtUtils::F77->runtime;

my @obj = @src;
map {s/\.[fc]/\.o/;} @obj; # swap .f, .c for .o

WriteMakefile(
	      'NAME'  	      => 'PDL::McIDAS',
	      'VERSION_FROM'  => 'McIDAS.pp',
	      

	      'TYPEMAPS'      => [&PDL_TYPEMAP()], 
	      'OBJECT'        => 'McIDAS.o ' . join (" ", @obj),
	      'DEFINE'        => $define,
	      'PM'	      => { 'McIDAS.pm' => '$(INST_LIBDIR)/McIDAS.pm'},
	      'INC'           => &PDL_INCLUDE(), # add include dirs as required by your lib
	      'LIBS'          => "$runtime -lg2c",  # add link directives as necessary
	      'clean'         => {'FILES'  => 'McIDAS.c McIDAS.pm McIDAS.xs McIDAS.o ' . join (" ", @obj) },
	      'dist'          => { COMPRESS => 'gzip', SUFFIX => 'gz' },       	
	      );

# Add genpp rule; this will invoke PDL::PP on our PP file
# the argument is an array reference where the array has three string elements:
#   arg1: name of the source file that contains the PP code
#   arg2: basename of the xs and pm files to be generated  
#   arg3: name of the package that is to be generated
sub MY::postamble { 
  $myfflags = ExtUtils::F77->cflags();
  my $orig = pdlpp_postamble(["McIDAS.pp",McIDAS,PDL::McIDAS]); 

  my $added = '';
  foreach my $s (@src) {
    my $o = $s;
    $o =~ s/\.[cf]/\.o/;
    
    if ($s =~ /\.c$/) { # C code
      $added .= "$o: $s\n\t$Config{'cc'} -c -o $o -g $s\n\n";
    } else { # FORTRAN code
      $added .= "$o: $s\n\t$myf77 -c -o $o -g $myfflags $underscore $s\n\n";
    }
    
  }

  return $orig . $added;

}  
