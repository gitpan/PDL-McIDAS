#include <stdio.h>
#include <stdlib.h>
#include <netinet/in.h>

#define PARAMETERS 9
#define GRIDS      3000


void l_grdlist(long int header[PARAMETERS][GRIDS], char file[1][256]) { 



  FILE* LECTURA;


  char label[33];
  long int project, date;
  long int max_grid;
  long int offset[GRIDS];

  long int elements, rows, columns;
  long int fecha, year, day, hour, validity;
  long int gridscale;
  long int gridname,gridunits,levelunits,maptype;
  long int level, levelscale, gridtype;
  long int timevariation, levelvariation, gridorigin;
  long int proj_34, proj_35, proj_36, proj_37;
  long int proj_38, proj_39;

  int n;

  int dato;
   
  printf("hola\n");
  printf(">%s<\n",kk[0]);
  printf("header %ld\n",*header);

    printf("header %d\n",header);


  dato=0;
  strcpy(archivo[0],"GRID0101");

  if ((LECTURA = fopen(archivo[0], "r")) == NULL) {
    printf("Error reading McIDAS file\n");
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
    (dato)++;
    fseek(LECTURA, offset[n]*4, SEEK_SET); 
    
    fread(&elements,4,1,LECTURA);
    elements=htonl(elements);
    
    fread(&rows,4,1,LECTURA);
    rows=htonl(rows);
    
    fread(&columns,4,1,LECTURA);
    columns=htonl(columns);
    
    fread(&fecha,4,1,LECTURA);
    fecha=htonl(fecha); 
    day=fecha%1000;
    year=((fecha+1900000)-day)/1000;     
    
    fread(&hour,4,1,LECTURA);
    hour=htonl(hour);
    fread(&validity,4,1,LECTURA);
    validity=htonl(validity);
    
    fread(&gridname,4,1,LECTURA);
    gridname=htonl(gridname);
    fread(&gridscale,4,1,LECTURA);
    gridscale=htonl(gridscale);
    fread(&gridunits,4,1,LECTURA);
    gridunits=htonl(gridunits);
    
    fread(&level,4,1,LECTURA);
    level=htonl(level);
    fread(&levelscale,4,1,LECTURA);
    levelscale=htonl(levelscale);
    fread(&levelunits,4,1,LECTURA);
    levelunits=htonl(levelunits);
    
    fread(&gridtype,4,1,LECTURA);
    gridtype=htonl(gridtype);
    fread(&timevariation,4,1,LECTURA);
    timevariation=htonl(timevariation);
    fread(&levelvariation,4,1,LECTURA);
    levelvariation=htonl(levelvariation);
    
    /* Variables reservadas */
    fseek(LECTURA, 16*4, SEEK_CUR);
    fread(&gridorigin,4,1,LECTURA);
    gridorigin=htonl(gridorigin);
    fread(&maptype,4,1,LECTURA);
    maptype=htonl(maptype);
    
    /* Variables de proyeccion */
    fread(&proj_34,4,1,LECTURA);
    proj_34=htonl(proj_34);
    fread(&proj_35,4,1,LECTURA);
    proj_35=htonl(proj_35);
    fread(&proj_36,4,1,LECTURA);
    proj_36=htonl(proj_36);
    fread(&proj_37,4,1,LECTURA);
    proj_37=htonl(proj_37);
    fread(&proj_38,4,1,LECTURA);
    proj_38=htonl(proj_38);
    fread(&proj_39,4,1,LECTURA);
    proj_39=htonl(proj_39);
    
    
    header[0][n] = gridname;
    header[1][n] = level;
    header[2][n] = levelunits;
    header[3][n] = year;
    header[4][n] = day;
    header[5][n] = hour;
    header[6][n] = validity;
    header[7][n] = (long int) n+1;
    header[8][n] = maptype;
  }
  
  fclose(LECTURA);   
}

for (n = 0; n < 5; n++) { 
  printf("n %d  ",n);
  printf("%4ld  %4ld %2ld ", header[0][n], header[1][n], header[2][n]);
  printf("%4ld %3ld %6ld  ", header[3][n], header[4][n], header[5][n]);
  printf("%3ld  %4ld %4ld\n", header[6][n], header[7][n], header[8][n]);  
}

printf("adios\n");


}




