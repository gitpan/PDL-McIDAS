# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.pl'

######################### We start with some black magic to print on failure.

# Change 1..1 below to 1..last_test_to_print .
# (It may become useful if the test is moved to ./t subdirectory.)

BEGIN { $| = 1; print "1..3\n"; }
END {print "not ok 1\n" unless $loaded;}
use PDL;
#use PDL::Char;
use PGPLOT;
use PDL::McIDAS;
use PDL::Graphics::PGPLOT;

$loaded=1;

######################### End of black magic.

# Insert your test code below (better if it prints "ok 13"
# (correspondingly "not ok 13") depending on the success of chunk 13
# of the test code):

#
# Convert date, lat, lon to sun fixed lat/lon
#


#$file=PDL::Char->new(zeroes(byte,256,1));
#$file->setstr(0,'GRID0101');
#$gridname=PDL::Char->new(zeroes(byte,4,3000));
#$levelunits=PDL::Char->new(zeroes(byte,4,3000));
#$maptype=PDL::Char->new(zeroes(byte,4,3000));
##@dfile=dims($file);
##print "dims file @dfile\n";
#
#$level=zeroes(long,3000);
#$rows=zeroes(long,3000);
#$columns=zeroes(long,3000);
#$year=zeroes(long,3000);
#$day=zeroes(long,3000);
#$hour=zeroes(long,3000);
#$validity=zeroes(long,3000);
#$num_grid=zeroes(long,3000);
#$num_grids=pdl(long,0);
#
#print "file    ",$file,"\n";
#
#h_grdlist($level,$rows,$columns,$year,$day,$hour,$validity,
#	   $num_grid,$num_grids,$file,$gridname,
#	   $levelunits,$maptype);
#
#print $level->slice('0:10'),"\n";
#print $rows->slice('0:10'),"\n";
#print $columns->slice('0:10'),"\n";
#print $year->slice('0:10'),"\n";
#print $hour->slice('0:10'),"\n";
#print $validity->slice('0:10'),"\n";
#print $num_grid->slice('0:10'),"\n";
#print $num_grids,"\n";
#print $file,"\n";
#print $gridname->slice(':,1:10'),"\n";
#print $levelunits->slice(':,1:10'),"\n";
#print $maptype->slice(':,1:10'),"\n";

    

($level,$rows,$columns,$year,$day,$hour,$validity,
 $num_grid,$maptype,$proj_34,$proj_35,$proj_36,
 $proj_37,$proj_38,$proj_39,
 $num_grids,$pfile,$gridname,
 $levelunits,$gridorigin)=PDL::McIDAS::grdlist('GRID0101');

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

print "ok 1\n";

#$file=PDL::Char->new(zeroes(byte,256,1));
#$file->setstr(0,'GRID0101');
#
#$data=zeroes(long,8456);
#
#$num_grid=pdl(long,5);
#
#h_grd2pdl($data,$num_grid,$file);
#
#print $data->slice('0:10');

$data=PDL::McIDAS::grd2pdl('GRID0101',5);

print $data->slice('0:10,0:10');

print "ok 2\n";

$ir=pdl(byte,[[255,255],[255,255]]);
$ig=pdl(byte,[[0,255],[255,0]]);
$ib=pdl(byte,[[0,0],[255,255]]);

dev('/xserve');
env(-0.5,1.5,-0.5,1.5);
$tr=pdl( [ 0, 1, 0, 0, 0, 1 ]);
env(0,1,0,1);
PDL::McIDAS::mcimagrgb($ir,$ig,$ib,{TRANSFORM=>$tr});

print "ok 3\n";
