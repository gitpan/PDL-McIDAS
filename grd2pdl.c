#include <stdio.h>
#include <netinet/in.h>

long int read_word(FILE *fich)
{ long int word;

  fread(&word,4,1,fich);
  word=htonl(word);
  return word;
}

void grd2pdl(long int elements, char archivo[1][256],
	     long int * ngrid, long int * data)
{ FILE* LECTURA;

  char label[33];
  long int project, date;
  long int max_grid;
  long int offset[3000];

  long int elementos; 
  long int filas, columnas;


  int i;
  int n;


  if ((LECTURA = fopen(archivo[0], "r")) != NULL)
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
    for(n=0;n<max_grid;n++) 
      offset[n]=htonl(offset[n]);

    /* Encabezamiento de cada GRID;
       Se mostraran los parametros clave del Grid indicado por pantalla*/

    if (offset[*ngrid-1]!=-1)
    { fseek(LECTURA, offset[*ngrid-1]*4, SEEK_SET); 

      elementos=read_word(LECTURA);
      filas=read_word(LECTURA);
      columnas=read_word(LECTURA);

      /* Variables reservadas */
      fseek(LECTURA, 61*4, SEEK_CUR);


      /* Asignacion de los valores a las variables de salida */
      //nrows=filas;
      //ncolumns=columnas;
      for (i=0;i<elements;i++)
          data[i]=read_word(LECTURA);
    }

    fclose(LECTURA);   
  }


}






