
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

  ip_mat *a, *b, *c, *d;
  Bitmap *photo;

  photo = bm_load( "flower.bmp" );
  a = bitmap_to_ip_mat( photo );
  bm_free ( photo );
/*  ----------------------------------------
*   PARTE 1:  Strutture dati, operazioni matematiche e gestione memoria
*   test delle funzioni concat, subset, copy e random
*   free e create saranno utilizzate sempre quindi non serve testarle
*   ----------------------------------------  */
/*
  printf("Test ip_mat_concat su dimensioni h e w\n");

  b = ip_mat_concat ( a, a, 0 );
  photo = ip_mat_to_bitmap( b );
  bm_save( photo, "testConcath.bmp" );
  bm_free( photo );
  ip_mat_free( b );

  b = ip_mat_concat ( a, a, 1 );
  photo = ip_mat_to_bitmap( b );
  bm_save( photo, "testConcatw.bmp" );
  bm_free( photo );
  ip_mat_free( b );

  printf("Le foto sono salvate come testConcath e testConcathw\n");
  printf("Test ip_mat_subset\n ");

  b = ip_mat_subset( a, 25, 125, 25, 125 );
  photo = ip_mat_to_bitmap( b );
  bm_save( photo, "testSubset.bmp" );
  bm_free( photo );
  ip_mat_free( b );

  printf("La foto è salvata come testSubset\n" );
  printf("Test ip_mat_copy\n ");

  b = ip_mat_copy( a );
  c = ip_mat_concat ( a, b, 1 );
  photo = ip_mat_to_bitmap( c );
  bm_save( photo, "testCopy.bmp" );
  bm_free( photo );
  ip_mat_free( b );
  ip_mat_free( c );

  printf("La foto è salvata come testCopy\n");

  printf("Test ip_mat_init_random\n");

  b = ip_mat_copy( a );
  ip_mat_init_random( b, (b->stat[0].mean+b->stat[1].mean+b->stat[2].mean)/3.0, 3.7 );
  photo = ip_mat_to_bitmap( b );
  bm_save( photo, "testInitRandom.bmp" );
  bm_free( photo );
  ip_mat_free( b );

  printf("La foto è salvata come testInitRandom\n");
*/
/*  ----------------------------------------
*   PARTE 2: Operazioni semplici su immagini
*   test delle funzioni brighten, gray_scale, corrupt e blend
*   ----------------------------------------  */
/*
  printf("Test ip_mat_brighten\n la foto è salvata come testBrighten\n");

  b = ip_mat_brighten( a, 75 );
  c = ip_mat_concat( a, b, 1 );
  photo = ip_mat_to_bitmap( c );
  bm_save( photo, "testBrighten.bmp" );
  bm_free( photo );
  ip_mat_free( b );
  ip_mat_free( c );

  printf("Test ip_mat_to_gray_scale\n la foto è salvata come testGrayScale\n");

  b = ip_mat_to_gray_scale( a );
  c = ip_mat_concat( a, b, 1 );
  photo = ip_mat_to_bitmap( c );
  bm_save( photo, "testGrayScale.bmp" );
  bm_free( photo );
  ip_mat_free( b );
  ip_mat_free( c );

  printf("Test ip_mat_corrupt\n la foto è salvata come testCorrupt\n");

  b = ip_mat_corrupt( a, 75 );
  clamp(b,0.0,255.0);
  c = ip_mat_concat( a, b, 1 );
  photo = ip_mat_to_bitmap( c );
  bm_save( photo, "testCorrupt.bmp" );
  bm_free( photo );
  ip_mat_free( b );
  ip_mat_free( c );

  printf("Test ip_mat_blend\n la foto è salvata come testBlend\n");

  photo = bm_load( "fullmoon.bmp" );
  b = bitmap_to_ip_mat( photo );
  bm_free( photo );
  c = ip_mat_subset ( b, 0, a->h, 0, a->w );
  ip_mat_free( b );
  b = ip_mat_blend( a, c, 0.5 );
  ip_mat_free( c );
  c = ip_mat_concat( a, b, 1 );
  photo = ip_mat_to_bitmap( c );
  bm_save( photo, "testBlend.bmp" );
  bm_free( photo );
  ip_mat_free( b );
  ip_mat_free( c );
*/
  /*  ----------------------------------------
  *   PARTE 3: Convoluzione e filtraggio
  *   test delle funzioni convolve, padding, sharpen, edge, emboss, average e gaussian
  *   le funzioni convolve e padding vengono usate dalle altre Funzioni
  *   le funzioni rescale e clamp non vengono testate perchè usate da funzioni nella parte 2
  *   ----------------------------------------  */

  printf("Test create_sharpen_filter\n");

  b = ip_mat_padding( a, 1, 1 );
  c = create_sharpen_filter();
  d = ip_mat_convolve ( b, c );
  clamp(d,0.0,255.0);
  ip_mat_free( b );
  b = ip_mat_concat( a, d, 1 );
  photo = ip_mat_to_bitmap( b );
  ip_mat_free( d );
  bm_save( photo, "testSharpen.bmp" );
  bm_free( photo );
  ip_mat_free( b );
  ip_mat_free( c );
  printf("la foto è salvata come testSharpen\n");
/*
  printf("Test create_edge_filter\n la foto è salvata come testEdge\n");

  b = ip_mat_padding( a, 1, 1 );
  c = create_edge_filter();
  d = ip_mat_convolve ( b, c );
  clamp(d,0.0,255.0);
  ip_mat_free( b );
  b = ip_mat_concat( a, d, 1 );
  ip_mat_free( d );
  photo = ip_mat_to_bitmap( b );
  bm_save( photo, "testEdge.bmp" );
  bm_free( photo );
  ip_mat_free( b );
  ip_mat_free( c );

  printf("Test create_emboss_filter\n la foto è salvata come testEmboss\n");

  b = ip_mat_padding( a, 1, 1 );
  c = create_emboss_filter();
  d = ip_mat_convolve ( b, c );
  clamp(d,0.0,255.0);
  ip_mat_free( b );
  b = ip_mat_concat( a, d, 1 );
  ip_mat_free( d );
  photo = ip_mat_to_bitmap( b );
  bm_save( photo, "testEmboss.bmp" );
  bm_free( photo );
  ip_mat_free( b );
  ip_mat_free( c );
*/
/*
  printf("Test create_average_filter\n la foto è salvata come testAverage\n");
  printf("dove arrivo?\n");
  b = create_average_filter( 3, 3, 1 );
  printf("qui1\n");
  c = ip_mat_padding( a, ((b->w*b->h)-1)/4, ((b->w*b->h)-1)/4 );
  printf("qui2\n");
  d = ip_mat_convolve ( a, b );
  printf("qui3\n");
  clamp(d,0.0,255.0);
  ip_mat_free( b );
  b = ip_mat_concat( a, d, 1 );
  printf("qui4\n");
  photo = ip_mat_to_bitmap( b );
  ip_mat_free( d );
  printf("qui5\n");
  bm_save( photo, "testAverage.bmp" );
  bm_free( photo );
  ip_mat_free( b );
  ip_mat_free( c );
  ip_mat_free ( a );
*/
  printf("Test create_gaussian_filter\n la foto è salvata come testGaussian\n");

  b = create_gaussian_filter( b->w, b->h, 1, 1 );
  b = ip_mat_padding( a, 1, 1 );
  d = ip_mat_convolve ( b, c );
  clamp(d,0.0,255.0);
  ip_mat_free( b );
  b = ip_mat_concat( a, d, 1 );
  ip_mat_free( d );
  photo = ip_mat_to_bitmap( b );
  bm_save( photo, "testGaussian.bmp" );
  bm_free( photo );
  ip_mat_free( b );
  ip_mat_free( c );
  ip_mat_free ( a );

  return 0;
}

/*  ----------------------------------------
*   inserire testo
*   ----------------------------------------  */
