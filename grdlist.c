#include <stdio.h>
#include <netinet/in.h>

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
	      long int maptype[GRIDS],
	      long int proj_34[GRIDS],
	      long int proj_35[GRIDS],
	      long int proj_36[GRIDS],
	      long int proj_37[GRIDS],
	      long int proj_38[GRIDS],
	      long int proj_39[GRIDS],
	      long int * p_num_grids,
	      char file[1][256], 
	      char gridname[GRIDS][4],
	      char levelunits[GRIDS][4],
	      char gridorigin[GRIDS][4] ) {



  FILE* LECTURA;


  char label[33];
  long int project, date;
  long int max_grid;
  long int offset[GRIDS];

  long int elements;
  long int fecha;
  long int gridscale;
  char gridunits[4];
  long int levelscale, gridtype;
  long int timevariation, levelvariation;

  int n;

  long int num_grids;
   
  num_grids=0;

  if ((LECTURA = fopen(file[0], "r")) == NULL) {
    printf("Error opening McIDAS file\n");
    return(1);
  }
  /* Offset de encabezamiento de archivo GRID;
     Determinacion del número de busqueda del primer GRID */

  fread(label,32,1,LECTURA);
  label[32]=0;
  
  fread(&project,4,1,LECTURA);
  project=htonl(project);
  fread(&date,4,1,LECTURA);
  date=htonl(date);
  
  /* Number maximum of grids */
  fread(&max_grid,4,1,LECTURA);
  max_grid=htonl(max_grid);
  max_grid=-max_grid;
  
  fread(offset,4,max_grid+1,LECTURA);
  for (n=0;n<=max_grid;n++) 
    offset[n]=htonl(offset[n]);
  
  
  /* Header for each GRID */
  for (n=0;n<(int)max_grid;n++) { 
    
    if (offset[n]==-1) continue;
    num_grids++;
    fseek(LECTURA, offset[n]*4, SEEK_SET); 
    
    fread(&elements,4,1,LECTURA);
    elements=htonl(elements);
    
    fread(&rows[n],4,1,LECTURA);
    rows[n]=htonl(rows[n]);
    
    fread(&columns[n],4,1,LECTURA);
    columns[n]=htonl(columns[n]);
    
    fread(&fecha,4,1,LECTURA);
    fecha=htonl(fecha); 
    day[n]=fecha%1000;
    year[n]=((fecha+1900000)-day[n])/1000;     
    
    fread(&hour[n],4,1,LECTURA);
    hour[n]=htonl(hour[n]);
    fread(&validity[n],4,1,LECTURA);
    validity[n]=htonl(validity[n]);
    
    fread(gridname[n],4,1,LECTURA);

    fread(&gridscale,4,1,LECTURA);
    gridscale=htonl(gridscale);

    fread(gridunits,4,1,LECTURA);
    
    fread(&level[n],4,1,LECTURA);
    level[n]=htonl(level[n]);

    fread(&levelscale,4,1,LECTURA);
    levelscale=htonl(levelscale);

    fread(levelunits[n],4,1,LECTURA);
    
    fread(&gridtype,4,1,LECTURA);
    gridtype=htonl(gridtype);
    fread(&timevariation,4,1,LECTURA);
    timevariation=htonl(timevariation);
    fread(&levelvariation,4,1,LECTURA);
    levelvariation=htonl(levelvariation);
    
    /* Variables reservadas */
    fseek(LECTURA, 17*4, SEEK_CUR);

    fread(gridorigin[n],4,1,LECTURA);

    fread(&maptype[n],4,1,LECTURA);
    maptype[n]=htonl(maptype[n]);
    
    /* Variables de proyeccion */
    fread(&proj_34[n],4,1,LECTURA);
    proj_34[n]=htonl(proj_34[n]);
    fread(&proj_35[n],4,1,LECTURA);
    proj_35[n]=htonl(proj_35[n]);
    fread(&proj_36[n],4,1,LECTURA);
    proj_36[n]=htonl(proj_36[n]);
    fread(&proj_37[n],4,1,LECTURA);
    proj_37[n]=htonl(proj_37[n]);
    fread(&proj_38[n],4,1,LECTURA);
    proj_38[n]=htonl(proj_38[n]);
    fread(&proj_39[n],4,1,LECTURA);
    proj_39[n]=htonl(proj_39[n]);
    

    num_grid[n]=(long int) n+1;
    
  }
  
  fclose(LECTURA);   

  *p_num_grids=num_grids;

  return(0);
  
//  for (n = 0; n < 5; n++) { 
//    printf("n %d  ",n);
//    printf("%4ld  %4ld %2ld ", header[0][n], header[1][n], header[2][n]);
//    printf("%4ld %3ld %6ld  ", header[3][n], header[4][n], header[5][n]);
//    printf("%3ld  %4ld %4ld\n", header[6][n], header[7][n], header[8][n]);  
//  }
//  
//  printf("adios\n");
  
  
}





