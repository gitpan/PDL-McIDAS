$VERSION = '0.12';

#-----------------------------------------------------------------------------------

pp_def	('sflatlon',
	Pars => 'date(); geolat(); geolon(); [o]sflat(); [o]sflon();',
	GenericTypes => [D],
	Code => 'coords($P(date),$P(geolat), $P(geolon), $P(sflat), $P(sflon));', 
	Doc  => <<'EOD');
Computes the sun-fixed lat and lon (in degrees), given a date/time
and geodetic lat/lon.

Inputs:
	date -- Date/time in 14 digit empress format (double)
             ex:  19970202230000 is Feb 2, 1997 at 11:00pm GMT
	lat, lon -- geodetic lat, lon in degrees (-90 to 90, -180 to 180)

Outputs:
	sflat, sflon -- sun-fixed latitude and longitude in degrees
EOD

#-----------------------------------------------------------------------------------

pp_def	('greatcirc',
	Pars => 'lat1(); lon1(); lat2(); lon2(); [o]dist();',
	GenericTypes => [D],
	Code => 'greatcirc($P(lat1),$P(lon1), $P(lat2), $P(lon2), $P(dist));', 
	Doc  => <<'EOD');
Computes the great circle distance in km given two lat/lon pairs in degrees.

Input and output variables are PDLs.  They can be scalars or N-dimensional arrays.
It is also possible to pass in lat1/lon1 as vectors and lat2/lon2 as scalar PDLs.
This will compute the distance between each location in lat1/lon1 and the lat2/lon2
scalar, returning the results in a vector PDL.
EOD

#-----------------------------------------------------------------------------------

pp_def	('trapezoid',
	Pars => 'x(n); y(n); [o]integral();',
	GenericTypes => [D],
	Code => 'int n_size; n_size = $SIZE(n); trapezoid(&n_size, $P(x), $P(y), $P(integral));',
	Doc  => <<'EOD');
Computes the integral of y=f(x) using the trapeozoid method.  X must
be monotonically increasing.
EOD

#-----------------------------------------------------------------------------------

pp_def	('lreg',
	Pars => 'x(n); y(n); [o]lreg(n); [o]c0(); [o]c1()',
	GenericTypes => [D],
	Code => 
	'int n_size; n_size = $SIZE(n); lreg(&n_size, $P(x), $P(y), $P(lreg), $P(c0), $P(c1));',
	Doc  => <<'EOD');	
Linear Regression (fitting points to a line).
($lreg, $c0, $c1) = Dutil::lreg ($x, $y);
$x -- independent variable
$y -- dependent variable
$lreg -- linear fit of y=f(x):  lreg=linear(x)
$c0, $c1 -- Linear fit parameters:  z = c0 + c1 * x
EOD

#-----------------------------------------------------------------------------------

pp_def	('lregslope',
	Pars => 'x(n); y(n); [o]lreg(n); [o]c1()',
	GenericTypes => [D],
	Code => 
	'int n_size; n_size = $SIZE(n); lregslope(&n_size, $P(x), $P(y), $P(lreg), $P(c1));',
	Doc  => <<'EOD');
Linear Regression (fitting points to a line).
This version finds only the slope of the line.  The y offset (c0) is
contrained to be zero, ie the line must pass through the origin.
($lreg, $c1) = Dutil::lregslope ($x, $y);
$x -- independent variable
$y -- dependent variable
$lreg -- linear fit of y=f(x):  lreg=linear(x)
$c1 -- Linear fit parameter:  z = c1 * x
EOD

#-----------------------------------------------------------------------------------

pp_def	('find_geop',
	Pars => 'zsfc(); temp(n); pres(n); e(n); [o]geop(n)',
	GenericTypes => [D],
	Code => 
	'int n_size; n_size = $SIZE(n); find_geop(&n_size, $P(zsfc), $P(temp), $P(pres), $P(e), $P(geop));',
	Doc  => <<'EOD');	
Compute a set of geopotential heights (in km) for a data set.
 Inputs:     1) alt -- Base geopotential height (km)
             2) t   -- Array of temperature values (C)
             3) p   -- Array of pressure values (mb)
             4) e   -- Array of water Vp values (mb, -999 if undef)
            
my $geop = Dutil::find_geop ($zsfc, $t, $p, $e);
EOD

#-----------------------------------------------------------------------------------

pp_def	('geo2msl',
	Pars => 'geop(); lat(); lon(); [o]geom()',
	GenericTypes => [D],
	Code => 
	'geo2msl($P(geop), $P(lat), $P(lon), $P(geom));',
	Doc  => <<'EOD');
Convert a list of geopotential altitudes to mean sea level altitudes.
 Inputs:     1) h   -- list of geopotential altitudes (in km);
             2) lat -- latitude  of profile in degrees.
             3) lon -- longitude of profile in degrees.

my $msl = Dutil::geo2msl ($geop, $lat, $lon);
EOD

#-----------------------------------------------------------------------------------

pp_def	('msl2geo',
	Pars => 'geom(); lat(); lon(); [o]geop()',
	GenericTypes => [D],
	Code => 
	'msl2geo($P(geom), $P(lat), $P(lon), $P(geop));',
	Doc  => <<'EOD');
Convert a list of geometric mean sea level altitudes to geopotential height
Return a list of geopotential altitudes, in km.

 Inputs:     1) z   -- list of geometric altitudes (in km);
             2) lat -- latitude  of profile in degrees.
             3) lon -- longitude of profile in degrees.

my $geop = Dutil::msl2geo ($msl, $lat, $lon);
EOD

#-----------------------------------------------------------------------------------

pp_def	('oned',
	Pars => 'x(); line(n); [o]interp()',
	GenericTypes => [D],
	Code => 
	'int n_size; n_size = $SIZE(n); oned(&n_size, $P(x), $P(line), $P(interp));'
);

#-----------------------------------------------------------------------------------	

pp_def	('dplog1',
	Pars => 'x(n); y(n); in(); [o]out();',
	GenericTypes => [D],
	Code => 'int n_size;
	         n_size = $SIZE(n); 
	         dplog1($P(x), $P(y), $in(), n_size, $P(out));',
	Doc  => <<'EOD');
Log-linear interpolation.  Returns $yint, corresponding to $xint as 
defined by $y = f($x).  Best for finding altitude if you know pressure.
$yint = Dutil::dplog1 ($x, $y, $xint);
$x -- independent variable (pressure)
$y -- dependent variable (altitude)
$xint -- point to interpolate from
$yint -- result of interpolation
EOD

#-----------------------------------------------------------------------------------	

pp_def	('dplinear',
	Pars => 'x(n); y(n); in(); [o]out();',
	GenericTypes => [D],
	Code => 'int n_size;
	         n_size = $SIZE(n); 
	         dplinear($P(x), $P(y), $in(), n_size, $P(out));',
	Doc  => <<'EOD');	
Linear interpolation.  Returns $yint, corresponding to $xint as 
defined by $y = f($x).  
$yint = Dutil::dplinear ($x, $y, $xint);
$x -- independent variable 
$y -- dependent variable 
$xint -- point to interpolate from
$yint -- result of interpolation
EOD

#-----------------------------------------------------------------------------------	

pp_def	('dpspline',
	Pars => 'x(n); y(n); in(nin); ylimit(); int m(); [o]out(nin);',
	GenericTypes => [D],
	Code => 'int n_size, nin_size;
	         n_size = $SIZE(n);
	         nin_size = $SIZE(nin); 
	         dpspline($P(x), $P(y), $P(in), n_size, nin_size, $ylimit(), $m(), $P(out));',
	Doc  => <<'EOD');	
Cubic spline interpolation.  Returns $yint, corresponding to $xint as 
defined by $y = f($x).  
$yint = Dutil::dpspline ($x, $y, $xint, $ylimit, $m);
$ylimit -- limit for second derivative.  
           Sometimes helps weed out sharp kinks
$m -- The number of points to base the spline on.  
      Can be set from 4 to the number of input points (size of x)
$x -- independent variable 
$y -- dependent variable 
$xint -- points to interpolate from
$yint -- results of interpolation
EOD

#-----------------------------------------------------------------------------------	

# This is a superfluous pp_def that I have left in as an example.  'greatcirc' can
# be made to do the same job that 'greatcirc1' was intended to do.
pp_def	('greatcirc1',
	Pars => 'lat(n); lon(n); lat2(); lon2(); [o]dist(n);',
	GenericTypes => [D],
	Code => '
	             int i, n_size;
	             n_size = $SIZE(n);
	             for (i=0;i<n_size;i++) {
   	               greatcirc(&$lat(n=>i),&$lon(n=>i), $P(lat2), $P(lon2), &$dist(n=>i));
		     }
	         '
);

pp_addpm({At => Top}, <<'EOD');
=head1 NAME

Dutil -- PDL interface to Doug Hunt's utility routines. Used at UCAR for
the GPS/MET and COSMIC projects.

=head1 SYNOPSIS

  use PDL;
  use Dutil;

  my ($date, $lat, $lon) = (1997020223000, 40.0, -105.0);

  my ($sflat, $sflon) = Dutil::sflatlon ($date, $lat, $lon);

  my ($lat1, $lon1 $lat2, $lon2) = (40.0, -105.0, 60.0, -130.0);

  my ($distance_km) = Dutil::greatcirc($lat1, $lon1, $lat2, $lon2);
  my ($integral)    = Dutil::trapezoid($x, $y); # where y = f(x)
  my ($lreg, $c0, $c1) = Dutil::lreg ($x, $y);	
  my ($lreg, $c1)   = Dutil::lregslope ($x, $y);	

=head1 DESCRIPTION 

required for pod2man

=head1 AUTHOR

	Doug Hunt (dhunt@ucar.edu)

=cut

EOD

pp_done();
