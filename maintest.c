
#include "ip_lib.h"
#include "bmp.h"
#include <stdlib.h>
#include <stdio.h>


int main(){
  ip_mat *a, *b, *c;
  Bitmap *photo;
  int i;
  photo = bm_load("flower.bmp");
  a = bitmap_to_ip_mat(photo);
  b = ip_mat_padding(a,1,1);
  /*bm_free( photo );
  photo = bm_load("vietnam.bmp");
  b = bitmap_to_ip_mat(photo);
*/

  /*c = ip_mat_blend( a,b,0.7);
  ip_mat_free( a );
  */
  ip_mat_free( a );
  a = ip_mat_create(3,3,1,1);
  for (i =0;i<3; i++)
  {
    set_val(a,i,1,0,0);
    set_val(a,i,2,0,-1);
  }
  ip_mat_show ( a );
  c = ip_mat_convolve(b,a);
  compute_stats( c );
  clamp (c,0,255);
  /*rescale(c,255);*/
  photo = ip_mat_to_bitmap( c );
  bm_save( photo, "testphoto.bmp" );
  bm_free( photo );

  ip_mat_free( a );
  ip_mat_free( b );

  return 0;
}
