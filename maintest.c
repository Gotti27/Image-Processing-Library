
#include "ip_lib.h"
#include "bmp.h"
#include <stdlib.h>
#include <stdio.h>


int main(){
  ip_mat *a, *b, *c;
  Bitmap *photo;

  photo = bm_load("flower.bmp");
  a = bitmap_to_ip_mat(photo);

  photo = bm_load("fullmoon.bmp");
  b = bitmap_to_ip_mat(photo);


  c = ip_mat_corrupt( a, 50.0);

  photo = ip_mat_to_bitmap( c );
  bm_save( photo, "testphoto.bmp" );
  bm_free( photo );

  ip_mat_free( a );

  return 0;
}

