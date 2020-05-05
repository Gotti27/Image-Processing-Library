#include "ip_lib.h"
#include "bmp.h"
#include <stdlib.h>
#include <stdio.h>


int main(){
  ip_mat *a;
  printf("\n\tStart your engines\n");
  a = ip_mat_create(5, 5, 3, 10.0);
  printf("Nope");


  ip_mat_show( a );

  /*printf("Hello World!\n");*/
  return 0;
}
