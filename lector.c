#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PARAMETERS 9
#define GRIDS      3000


int h_grdlist(long int level[GRIDS],
	      long int rows[GRIDS],
	      long int columns[GRIDS],
	      long int year[GRIDS],
	      long int day[GRIDS],
	      long int hour[GRIDS],
	      long int validity[GRIDS],
	      long int num_grid[GRIDS],
	      long int * p_num_grids,
	      char file[1][256],
	      char gridname[GRIDS][4],
	      char levelunits[GRIDS][4],
	      char maptype[GRIDS][4]);


int main(int argc, char* argv[]) { 

  long int level[GRIDS];
  long int rows[GRIDS];
  long int columns[GRIDS];
  long int year[GRIDS];
  long int day[GRIDS];
  long int hour[GRIDS];
  long int validity[GRIDS];
  long int num_grid[GRIDS];
  long int num_grids;
  char file[1][256];
  char gridname[GRIDS][4];
  char levelunits[GRIDS][4];
  char maptype[GRIDS][4];


  int i;

  strcpy(file[0],"GRID0101");



  h_grdlist(level,
	    rows,
	    columns,
	    year,
	    day,
	    hour,
	    validity,
	    num_grid,
	    &num_grids,
	    file,
	    gridname,
	    levelunits,
	    maptype);
  
  

  
  for (i = 0; i < 10; i++) { 
    //printf("i %d\n",i);
    printf("%ld  %ld %ld ", level[i], rows[i], columns[i]);
    printf("%ld %ld %ld  ", year[i], day[i], hour[i]);
    printf("%ld  %ld\n", validity[i], num_grid[i]);  
  }
  
  
  exit(0);
  
}





