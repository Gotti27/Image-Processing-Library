
#include "ip_lib.h"
#include "bmp.h"
#include <stdlib.h>
#include <stdio.h>

/*  ----------------------------------------
*   Test funzioni del progetto V0
*
*   Ricordarsi che :
*   per liberare una bitmap -> bm_free( bitmap );
*   per liberare una ip_mat -> ip_mat_free( ip_mat );
*
*   ----------------------------------------
*   Funzioni Libreria bmp.h utilizzabili
*
*   bm_create( larghezza , altezza );
*   bm_load( "immagine.bmp" );
*   bm_save( bitmap , "nomeFotoSalvata.bmp" );
*   bm_free( bitmap );
*
*   ----------------------------------------  */


int main(){

  ip_mat *a, *b, *c;
  Bitmap *photo;
  int i,j,l;

/*  ----------------------------------------
*   Parte 1:  Strutture dati, operazioni matematiche e gestione memoria
*   test delle funzioni concat, subset, copy e random
*   free e create saranno utilizzate sempre quindi non serve testarle
*   ----------------------------------------  */

  printf("Test ip_mat_concat su dimensioni h e w\n le foto sono salvate come testConcath e testConcathw\n");
  
  photo = bm_load( "flower.bmp" );
  a = bitmap_to_ip_mat( photo );
  
  b = ip_mat_concat ( a , a , 0 );
  photo = ip_mat_to_bitmap( b );
  bm_save( photo, "testConcath.bmp" );
  bm_free( photo );
  ip_mat_free( b );

  b = ip_mat_concat ( a , a , 1 );
  photo = ip_mat_to_bitmap( b );
  bm_save( photo, "testConcatw.bmp" );
  bm_free( photo );
  ip_mat_free( b );
  

  printf("Test ip_mat_subset\n la foto è salvata come testSubset\n");
  
  b = ip_mat_subset( a , 25 , 125 , 25 , 125 );
  photo = ip_mat_to_bitmap( b );
  bm_save( photo, "testSubset.bmp" );
  bm_free( photo );
  ip_mat_free( b );

  printf("Test ip_mat_copy\n la foto è salvata come testCopy\n");

  b = ip_mat_copy( a );
  c = ip_mat_concat ( a , b , 1);
  photo = ip_mat_to_bitmap( c );
  bm_save( photo, "testCopy.bmp" );
  bm_free( photo );
  ip_mat_free( b );
  ip_mat_free( c );

  printf("Test ip_mat_init_random\n la foto è salvata come testInitRandom\n");

  b = ip_mat_create( 100 , 100 , 3 , 0 );
  ip_mat_init_random( b , 100 , 3.7 );
  photo = ip_mat_to_bitmap( b );
  bm_save( photo, "testInitRandom.bmp" );
  bm_free( photo );
  ip_mat_free( b );

/*  ----------------------------------------
*   PARTE 2: Operazioni semplici su immagini
*   test delle funzioni brighten, gray_scale, corrupt e blend
*   ----------------------------------------  */

  printf("Test ip_mat_brighten\n la foto è salvata come testBrighten\n");

  b = ip_mat_brighten( a , 75);
  c = ip_mat_concat( a , b , 1);
  photo = ip_mat_to_bitmap( c );
  bm_save( photo, "testBrighten.bmp" );
  bm_free( photo );
  ip_mat_free( b );
  ip_mat_free( c );

  printf("Test ip_mat_gray_scale\n la foto è salvata come testGrayScale\n");

  b = ip_mat_gray_scale( a );
  c = ip_mat_concat( a , b , 1);
  photo = ip_mat_to_bitmap( c );
  bm_save( photo, "testGrayScale.bmp" );
  bm_free( photo );
  ip_mat_free( b );
  ip_mat_free( c ); 

  printf("Test ip_mat_corrupt\n la foto è salvata come testCorrupt\n");

  b = ip_mat_corrupt( a , 125 );
  c = ip_mat_concat( a , b , 1);
  photo = ip_mat_to_bitmap( c );
  bm_save( photo, "testCorrupt.bmp" );
  bm_free( photo );
  ip_mat_free( b );
  ip_mat_free( c ); 

  printf("Test ip_mat_blend\n la foto è salvata come testBlend\n");

  photo = bm_load( "moon.bmp" );
  b = bitmap_to_ip_mat( photo );
  bm_free( photo ):
  c = ip_mat_subset ( b, 0, a->w, 0, a->h );
  ip_mat_free( b );
  b = ip_mat_blend( a, c, 0.5 );
  ip_mat_free( c );
  c = ip_mat_concat( a , c , 1);
  photo = ip_mat_to_bitmap( c );
  bm_save( photo, "testBlend.bmp" );
  bm_free( photo );
  ip_mat_free( b );
  ip_mat_free( c ); 


  ip_mat_free ( a );

  return 0;
}

/*  ----------------------------------------
*   inserire testo  
*   ----------------------------------------  */