#!/usr/bin/perl -w

use PDL;
use PGPLOT;

@x=(1,2,3,4);
@y=(5,1,2,7);

pgbegin(0,"/xserve",1,1);
pgenv(1,5,1,10,0,0);
pgpoint(4,\@x,\@y,0);
