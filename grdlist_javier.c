#include <stdio.h>
#include <stdlib.h>
#include <netinet/in.h>

unsigned long int read_word(FILE *fich)
{ unsigned long int word;

  fread(&word,4,1,fich);
  word=htonl(word);
  return word;
}

main(int argc, char *argv[])
{ FILE *LECTURA;

  char label[33];
  unsigned long int project, date;
  long int max_grid;
  unsigned long int offset[3000];

  unsigned long int elementos, filas, columnas;
  unsigned long int fecha, anno, dia, hora, validez;
  unsigned long int gridscale;
  char gridname[5],gridunits[5],levelunits[5],maptype[5];
  unsigned long int level, levelscale, gridtype;
  unsigned long int timevariation, levelvariation, gridorigin;
  unsigned long int proj_34, proj_35, proj_36, proj_37;
  long int proj_38, proj_39;
  long int dato;

  int n;

  if ((argc==2) && ((LECTURA = fopen(argv[1], "r")) != NULL))
  {  /* Offset de encabezamiento de archivo GRID;
        Determinacion del número de busqueda del primer GRID */
  
    fread(label,32,1,LECTURA);
    label[32]=0;

    project=read_word(LECTURA);  
    date=read_word(LECTURA);

    fread(&max_grid,4,1,LECTURA);
    max_grid=htonl(max_grid);
    max_grid=-max_grid;
      
    fread(offset,4,max_grid+1,LECTURA);
    for (n=0;n<=max_grid;n++) 
      offset[n]=htonl(offset[n]);

    /* Encabezamiento de cada GRID;
       Se mostraran por pantalla los parametros clave de cada Grid */
//  printf("PARAM  LEVEL    DAY     TIME  FHOUR GRID PROJ\n");
//  printf("----- ------- -------- ------ ----- ---- ----\n");   

    for (n=0;n<max_grid;n++)
    { if (offset[n]==-1) continue;
      fseek(LECTURA, offset[n]*4, SEEK_SET); 

      elementos=read_word(LECTURA);
      filas=read_word(LECTURA);
      columnas=read_word(LECTURA);

      fecha=read_word(LECTURA); 
      dia=fecha%1000;
      anno=((fecha+1900000)-dia)/1000;     

      hora=read_word(LECTURA);
      validez=read_word(LECTURA);

      fread(gridname,4,1,LECTURA);	gridname[4]=0;
      gridscale=read_word(LECTURA);
      fread(gridunits,4,1,LECTURA);	gridunits[4]=0;

      level=read_word(LECTURA);
      levelscale=read_word(LECTURA);
      fread(levelunits,4,1,LECTURA);	levelunits[4]=0;
       
      gridtype=read_word(LECTURA);
      timevariation=read_word(LECTURA);
      levelvariation=read_word(LECTURA);

      /* Variables reservadas */
      fseek(LECTURA, 16*4, SEEK_CUR);
      gridorigin=read_word(LECTURA);
      fread(maptype,4,1,LECTURA);	maptype[4]=0;

      /* Variables de proyeccion */
      proj_34=read_word(LECTURA);
      proj_35=read_word(LECTURA);
      proj_36=read_word(LECTURA);
      proj_37=read_word(LECTURA);
      fread(&proj_38,4,1,LECTURA);	proj_38=htonl(proj_38);
      fread(&proj_39,4,1,LECTURA);	proj_39=htonl(proj_39);

      printf("%4s  %4u %c%c ", gridname, level, levelunits[0], levelunits[1]);
      printf("%4d %3d %6d  %3d  ", anno, dia, hora, validez); 
      printf("%4d %4s\n", n+1, maptype);       
    }

    fclose(LECTURA);   
  }
}





