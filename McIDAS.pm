
#
# GENERATED WITH PDL::PP! Don't modify!
#
package PDL::McIDAS;

@EXPORT_OK  = qw( PDL::PP h_grdlist PDL::PP h_grd2pdl );
%EXPORT_TAGS = (Func=>[@EXPORT_OK]);

use PDL::Core;
use PDL::Exporter;
use DynaLoader;



   
   @ISA    = ( 'PDL::Exporter','DynaLoader' );
   push @PDL::Core::PP, __PACKAGE__;
   bootstrap PDL::McIDAS ;




=head1 NAME

PDL::McIDAS -- PDL interface to McIDAS.

=head1 SYNOPSIS

  use PDL;
  use PDL::McIDAS;

  ($level,$rows,$columns,$year,$day,$hour,$validity,
 	$num_grid,$maptype,$proj_34,$proj_35,$proj_36,
	 $proj_37,$proj_38,$proj_39,
	 $num_grids,$pfile,$gridname,
	 $levelunits,$gridorigin)=
            PDL::McIDAS::grdlist('GRID0101');

  $data=PDL::McIDAS::grd2pdl('GRID0101',5);

  PDL::McIDAS::mcimagrgb($ir,$ig,$ib,{TRANSFORM=>$tr});

=head1 DESCRIPTION 

PDL::McIDAS allows you to read into PDL the McIDAS GRID files. There
is also a function to plot RGB images.

=head1 AUTHOR

	Xavier Calbet    and Javier Garcia-Pereda 
        (xcalbet@yahoo.es, javierpereda@yahoo.com)

Copyright (c) 2003 Xavier Calbet and Javier Garcia-Pereda 
(xcalbet@yahoo.es, javierpereda@yahoo.com).  All
rights reserved. This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

=head1 FUNCTIONS

=head2 grdlist

PDL::McIDAS::grdlist is equivalent to the GRDLIST command from McIDAS
and shows the various header parameters of the McIDAS GRID file.
It accepts as argument a McIDAS GRID file name. As output it gives an array
of PDLs, each one of them containing a different parameter of the GRID
header. As an example, if a McIDAS GRID file called "GRID0101" 
is to be read into PDL we will proceed as follows:

($level,$rows,$columns,$year,$day,$hour,$validity,
 $num_grid,$maptype,$proj_34,$proj_35,$proj_36,
 $proj_37,$proj_38,$proj_39,
 $num_grids,$pfile,$gridname,
 $levelunits,$gridorigin)=
        PDL::McIDAS::grdlist('GRID0101');

$cgridname=PDL::Char->new(PDL::copy($gridname));

print "level      ",$level,"\n";

print "rows       ",$rows,"\n";

print "columns    ",$columns,"\n";

print "year       ",$year,"\n";

print "day        ",$day,"\n";

print "hour       ",$hour,"\n";

print "validity   ",$validity,"\n";

print "num_grid   ",$num_grid,"\n";

print "maptype    ",$maptype,"\n";

print "proj_34    ",$proj_34,"\n";

print "proj_35    ",$proj_35,"\n";

print "proj_36    ",$proj_36,"\n";

print "proj_37    ",$proj_37,"\n";

print "proj_38    ",$proj_38,"\n";

print "proj_39    ",$proj_39,"\n";

print "num_grids  ",$num_grids,"\n";

print "pfile      ",$pfile,"\n";

print "gridname   ",$gridname,"\n";

print "cgridname  ",$cgridname,"\n";

print "levelunits ",$levelunits,"\n";

print "gridorigin ",$gridorigin,"\n";

=head2 grd2pdl

PDL::McIDAS::grd2pdl reads the data from a McIDAS grid and transfers
it to a PDL. As arguments it accepts a McIDAS GRID file name and
the grid number to be read. As output, a PDL is given containing
the data from the GRID file. The grid number can be found with the
grdlist function in the $num_grid PDL header parameter of the McIDAS
GRID file (see the description of the grdlist function in this man page).
As an example:

$data=PDL::McIDAS::grd2pdl('GRID0101',5);

=head2 mcimagrgb

PDL::McIDAS::mcimagrgb draws an RGB image with PGPLOT, composing
colors with the RGB data. As input it accepts three two dimensional
 PDLs with ranges between 0 and 255 which correspond to the red, green
and blue colors and the modifiers to the PGPLOT imag function.
As an example:

$ir=pdl(byte,[[255,255],[255,255]]);
$ig=pdl(byte,[[0,255],[255,0]]);
$ib=pdl(byte,[[0,0],[255,255]]);

dev('/xserve');

env(-0.5,1.5,-0.5,1.5);

$tr=pdl( [ 0, 1, 0, 0, 0, 1 ]);

env(0,1,0,1);

PDL::McIDAS::mcimagrgb($ir,$ig,$ib,{TRANSFORM=>$tr});

=cut


sub grdlist {
    
    use PDL::Char;
    
    ($file)=@_;
    
    $pfile=PDL::Char->new(zeroes(byte,256,1));
    $pfile->setstr(0,$file);
    $gridname=zeroes(byte,4,3000);
    $levelunits=zeroes(byte,4,3000);
    $gridorigin=zeroes(byte,4,3000);
    
    $level=zeroes(long,3000);
    $rows=zeroes(long,3000);
    $columns=zeroes(long,3000);
    $year=zeroes(long,3000);
    $day=zeroes(long,3000);
    $hour=zeroes(long,3000);
    $validity=zeroes(long,3000);
    $num_grid=zeroes(long,3000);
    $maptype=zeroes(long,3000);
    $proj_34=zeroes(long,3000);
    $proj_35=zeroes(long,3000);
    $proj_36=zeroes(long,3000);
    $proj_37=zeroes(long,3000);
    $proj_38=zeroes(long,3000);
    $proj_39=zeroes(long,3000);
    $num_grids=pdl(long,0);
    
    h_grdlist($level,$rows,$columns,$year,$day,$hour,$validity,
	      $num_grid,$maptype,$proj_34,$proj_35,$proj_36,
	      $proj_37,$proj_38,$proj_39,
	      $num_grids,$pfile,$gridname,
	      $levelunits,$gridorigin);


    # Reducing the dimensionality to save memory
    $ngridname=PDL::copy($gridname->slice(':,0:'.($num_grids-1)));
    $gridname=null;
    $nlevelunits=PDL::copy($levelunits->slice(':,0:'.($num_grids-1)));
    $levelunits=null;
    $ngridorigin=PDL::copy($gridorigin->slice(':,0:'.($num_grids-1)));
    $gridorigin=null;
    $nlevel=PDL::copy($level->slice('0:'.($num_grids-1)));
    $level=null;
    $nrows=PDL::copy($rows->slice('0:'.($num_grids-1)));
    $rows=null;
    $ncolumns=PDL::copy($columns->slice('0:'.($num_grids-1)));
    $columns=null;
    $nyear=PDL::copy($year->slice('0:'.($num_grids-1)));
    $year=null;
    $nday=PDL::copy($day->slice('0:'.($num_grids-1)));
    $day=null;
    $nhour=PDL::copy($hour->slice('0:'.($num_grids-1)));
    $hour=null;
    $nvalidity=PDL::copy($validity->slice('0:'.($num_grids-1)));
    $validity=null;
    $nnum_grid=PDL::copy($num_grid->slice('0:'.($num_grids-1)));
    $num_grid=null;
    $nmaptype=PDL::copy($maptype->slice('0:'.($num_grids-1)));
    $maptype=null;
    $nproj_34=PDL::copy($proj_34->slice('0:'.($num_grids-1)));
    $proj_34=null;
    $nproj_35=PDL::copy($proj_35->slice('0:'.($num_grids-1)));
    $proj_35=null;
    $nproj_36=PDL::copy($proj_36->slice('0:'.($num_grids-1)));
    $proj_36=null;
    $nproj_37=PDL::copy($proj_37->slice('0:'.($num_grids-1)));
    $proj_37=null;
    $nproj_38=PDL::copy($proj_38->slice('0:'.($num_grids-1)));
    $proj_38=null;
    $nproj_39=PDL::copy($proj_39->slice('0:'.($num_grids-1)));
    $proj_39=null;

    
    return ($nlevel,$nrows,$ncolumns,$nyear,$nday,$nhour,$nvalidity,
	    $nnum_grid,$nmaptype,$nproj_34,$nproj_35,$nproj_36,
	    $nproj_37,$nproj_38,$nproj_39,
	    $num_grids,$pfile,$ngridname,
	    $nlevelunits,$ngridorigin);
    
}


sub grd2pdl {
    
    use PDL::Char;
    
    ($file,$snum_grid)=@_;
    
    $pfile=PDL::Char->new(zeroes(byte,256,1));
    $pfile->setstr(0,$file);

    $pnum_grid=pdl(long,$snum_grid);

    ($level,$rows,$columns,$year,$day,$hour,$validity,
     $num_grid,$maptype,$proj_34,$proj_35,$proj_36,
     $proj_37,$proj_38,$proj_39,$num_grids,$pfile,$gridname,
     $levelunits,$gridorigin)=grdlist($file);

    $srows=at($rows,$snum_grid);
    $scols=at($columns,$snum_grid);

    $elements=$srows*$scols;

    $data=zeroes(long,$elements);

    h_grd2pdl($data,$pnum_grid,$pfile);
    
    $odata=$data->reshape($srows,$scols);

    return ($odata);
    
}

sub mcimagrgb {
    
    use strict;
    use PGPLOT;
    use PDL::ImageRGB;
    use PDL::Graphics::PGPLOT;


    my ($ir,$ig,$ib,$rh)=@_;

    my ($lo,$hi);
    my ($oldlo,$oldhi);
    my ($nlevelm1,$imag,$nimag);
    my ($out,$lut);
    my ($r,$g,$b);
    my ($levels);
    my ($i,@cr,@cg,@cb);

    # Remembering old colors to restore them
    #pgqcir($oldlo,$oldhi);
    #print "oldlo $oldlo $oldhi\n";
    #for ($i=0;$i<$oldlo;$i++) {
    #	pgqcr($i,$cr[$i],$cg[$i],$cb[$i]);
    #}
    #for ($i=0;$i<$oldlo;$i++) {
    #	print "color $i $cr[$i]  $cg[$i]  $cb[$i]\n";
    #}


    #pgqcol($lo,$hi);
    pgqcir($lo,$hi);
    #print "hi $hi lo $lo\n";
    if ($hi == 0 && $lo == 0) {
	barf("No PGPLOT device has been called\n");
    }
    
    $nlevelm1=($hi-$lo)*1.;
    $imag=cat($ir,$ig,$ib)->xchg(0,2)->xchg(1,2);
    $nimag=byte(
	      PDL::copy((($imag-min($imag))/
			 (max($imag)-min($imag))*($nlevelm1))));
    
    ($out,$lut)=cquant($nimag,$nlevelm1+1);
    $out=$out/$nlevelm1;
    #print "out $out\n";
    
    $lut=$lut/$nlevelm1;
    $r=$lut->slice('(0),:');
    $g=$lut->slice('(1),:');
    $b=$lut->slice('(2),:');
    
    $levels=sequence(float,($nlevelm1+1))/$nlevelm1;
    


    #pgscir($lo,$lo+$nlevelm1);
    ctab($levels,$r,$g,$b);

    imag($out,0,1,$rh);

    # Restoring old colors
    #for ($i=0;$i<$oldlo;$i++) {
    #	pgscr($i,$cr[$i],$cg[$i],$cb[$i]);
    #}
    #pgscir($oldlo,$oldhi);
}









=head1 FUNCTIONS



=cut






=head2 h_grdlist

=for sig

  Signature: ([o]level(3000);[o]rows(3000);[o]columns(3000);[o]year(3000);[o]day(3000);[o]hour(3000);[o]validity(3000);[o]num_grid(3000);[o]maptype(3000);[o]proj_34(3000);[o]proj_35(3000);[o]proj_36(3000);[o]proj_37(3000);[o]proj_38(3000);[o]proj_39(3000);[o]num_grids();byte file(256,1);byte [o]gridname(4,3000);byte [o]levelunits(4,3000);byte [o]gridorigin(4,3000))

Inner function (not oriented to the end user) 
used to lists the grid headers of the McIDAS grid file 
defined in grid_file.





=cut






*h_grdlist = \&PDL::h_grdlist;




=head2 h_grd2pdl

=for sig

  Signature: ([o]data(elements); num_grid(); byte file(256,1))

Inner function (not intended for end users) which copies a 
given grid from a McIDAS grid file to a pdl matrix.




=cut






*h_grd2pdl = \&PDL::h_grd2pdl;


;



# Exit with OK status

1;

		   