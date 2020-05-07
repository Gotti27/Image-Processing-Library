/*
 Created by Sebastiano Vascon on 23/03/20.
*/

#include <stdio.h>
#include "ip_lib.h"
#include "bmp.h"

void ip_mat_show(ip_mat * t){
    unsigned int i,l,j;
    printf("Matrix of size %d x %d x %d (hxwxk)\n",t->h,t->w,t->k);
    for (l = 0; l < t->k; l++) {
        printf("Slice %d\n", l);
        for(i=0;i<t->h;i++) {
            for (j = 0; j < t->w; j++) {
                printf("%f ", get_val(t,i,j,l));
            }
            printf("\n");
        }
        printf("\n");
    }
}

void ip_mat_show_stats(ip_mat * t){
    unsigned int k;

    compute_stats(t);

    for(k=0;k<t->k;k++){
        printf("Channel %d:\n", k);
        printf("\t Min: %f\n", t->stat[k].min);
        printf("\t Max: %f\n", t->stat[k].max);
        printf("\t Mean: %f\n", t->stat[k].mean);
    }
}

ip_mat * bitmap_to_ip_mat(Bitmap * img){
    unsigned int i=0,j=0;

    unsigned char R,G,B;

    unsigned int h = img->h;
    unsigned int w = img->w;

    ip_mat * out = ip_mat_create(h, w,3,0);

    for (i = 0; i < h; i++)              /* rows */
    {
        for (j = 0; j < w; j++)          /* columns */
        {
            bm_get_pixel(img, j,i,&R, &G, &B);
            set_val(out,i,j,0,(float) R);
            set_val(out,i,j,1,(float) G);
            set_val(out,i,j,2,(float) B);
        }
    }

    return out;
}

Bitmap * ip_mat_to_bitmap(ip_mat * t){

    Bitmap *b = bm_create(t->w,t->h);

    unsigned int i, j;
    for (i = 0; i < t->h; i++)              /* rows */
    {
        for (j = 0; j < t->w; j++)          /* columns */
        {
            bm_set_pixel(b, j,i, (unsigned char) get_val(t,i,j,0),
                    (unsigned char) get_val(t,i,j,1),
                    (unsigned char) get_val(t,i,j,2));
        }
    }
    return b;
}

float get_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k){
    if(i<a->h && j<a->w &&k<a->k){  /* j>=0 and k>=0 and i>=0 is non sense*/
        return a->data[i][j][k];
    }else{
        printf("Errore get_val!!!");
        exit(1);
    }
}

void set_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k, float v){
    if(i<a->h && j<a->w &&k<a->k){
        a->data[i][j][k]=v;
    }else{
        printf("Errore set_val!!!");
        exit(1);
    }
}

float get_normal_random(){
    float y1 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    float y2 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    return cos(2*PI*y2)*sqrt(-2.*log(y1));

}


ip_mat * ip_mat_create(unsigned int h, unsigned int w, unsigned int k, float v){
  unsigned int i,j,l;
  ip_mat *nuova;
  stats * st;
  float ***p3;

  printf("ok\n");

  nuova = malloc(sizeof(ip_mat));
  nuova->w = w;
  nuova->h = h;
  nuova->k = k;

  printf("ok\n");

  st = malloc(k * sizeof(stats));
  printf("Mallocato\n");
  for ( i = 0; i < k; i++){
    printf("%f", v);
    (&st[i])->min =  v;
    printf("okl");
    (&st[i])->max =  v;
    (&st[i])->mean = v;
    printf("okl\n");
  }
  printf("%f", v);
  nuova->stat = st;
  printf("ok\n");
  /* parkour */
  p3 = malloc (h * sizeof(float**));
  for ( i = 0; i < h; i++)
  {
    p3[i] =  malloc (w * sizeof(float*));
    for ( j = 0; j < w; j++)
    {
      p3[i][j] = malloc (k * sizeof(float));
      for ( l = 0; l < k; l++)
      {
        p3[i][j][l] = v;
      }
    }
  }
  printf("ok\n");
  nuova->data = p3;
  printf("ok\n");
  return nuova;
}

void compute_stats(ip_mat * t)
{
  unsigned int i,j,l;
  float min,max,mean;
  for(i=0;i< t->k;i++)
  {
    min =  t->data[0][0][i];/*penso qui scoppia*/
    max =  t->data[0][0][i];
    mean = 0.0;
  }
  for (i = 0; i < t->k; i++)
  {
    for ( j = 0; j < t->h; j++)
    {
      for ( l = 0; l < t->w; l++)
      {
        /*qui bisogna fare conti e ci sarà casino*/
        if(t->data[j][l][i] > max )
        {
          max = t->data[j][l][i];
        }
        else
        {
          if(t->data[j][l][i] < min )
          {
            min = t->data[j][l][i];
          }
        }
        mean += t->data[j][l][i] / (t->h * t->w);


      }
    }
    (&t->stat[i])->min = min;
    (&t->stat[i])->max = max;
    (&t->stat[i])->mean = mean;
  }
  printf("ho finito\n" );
}


void ip_mat_free(ip_mat * t){
  int i,j;
  free(t->stat);

  for (i = 0; i< t->h; i++){
    for (j = 0; j< t->w; j++){
      free(t->data[i][j]);
    }
    free(t->data[i]);
  }
  free(t->data);

  free(t);
}

void ip_mat_init_random(ip_mat * t, float mean, float var){
  int i, j, l;
  for (i = 0; i < t->h; i++){
    for (j = 0; j < t->w; j++){
      for (l = 0; l < t->k; l++){
        t->data[i][j][l] = var*get_normal_random() + mean;
      }
    }
  }
  compute_stats(t);
}

ip_mat * ip_mat_sum(ip_mat * a, ip_mat * b){
  if(a->w!=b->w || a->h!=b->h || a->k!=b->k){
      exit(1);
  }
  else{
    int i, j, l;
    ip_mat* sum = ip_mat_copy( a );

    for (i = 0; i < a->h; i++){
      for (j = 0; j < a->w; j++){
        for (l = 0; l < a->k; l++){
          sum->data[i][j][l] = a->data[i][j][l] + b->data[i][j][l];
        }
      }
    }
    compute_stats(sum);
    return sum;
  }

}

ip_mat * ip_mat_copy(ip_mat * in){
    int i,j,l;
    ip_mat * copia;
    copia = ip_mat_create(in->h, in->w, in->k, 1.0);
    for (i = 0; i < in->h; i++) {
        for(j = 0; j < in->w; j++){
            for(l = 0; l < in->k; l++){
                copia->data[i][j][l] = in->data[i][j][l];
            }
        }
    }
    compute_stats(copia);
    return copia;
}

ip_mat * ip_mat_sub(ip_mat * a, ip_mat * b){
    if (a->w == b->w && a->h == b->h && a->k == b->k) {
        int i,j,l;
        ip_mat * sub;
        sub = ip_mat_copy(a);
        for (i = 0; i < a->h; i++) {
            for(j = 0; j < a->w; j++){
                for(l = 0; l < a->k; l++){
                    sub->data[i][j][l] = a->data[i][j][l] - b->data[i][j][l];
                }
            }
        }
        compute_stats(sub);
        return sub;
    }
    else{
        exit(1);
    }
}

ip_mat * ip_mat_mean(ip_mat * a, ip_mat * b){
    if(a->w == b->w && a->h == b->h && a->k == b->k){
        int i,j,l;
        ip_mat * mean = ip_mat_copy( a );

        for (i = 0; i < a->h; i++) {
            for(j = 0; j < a->w; j++){
                for(l = 0; l < a->k; l++){
                    mean->data[i][j][l] = (a->data[i][j][l] + b->data[i][j][l]) / 2;
                }
            }
        }
        compute_stats(mean);
        return mean;
    }
    else{
        exit(1);
    }
}

ip_mat * ip_mat_sub(ip_mat * a, ip_mat * b){
    if (a->w == b->w && a->h == b->h && a->k == b->k) {
        int i,j,l;
        ip_mat * sub;
        sub = ip_mat_copy(a);
        for (i = 0; i < a->h; i++) {
            for(j = 0; j < a->w; j++){
                for(l = 0; l < a->k; l++){
                    sub->data[i][j][l] = a->data[i][j][l] - b->data[i][j][l];
                }
            }
        }
        compute_stats(sub);
        return sub;
    }
    else{
        exit(1);
    }
}

ip_mat * ip_mat_mul_scalar(ip_mat *a, float c){
 ip_mat mus = ip_mat_create(a->h,a->w,a->k,0);
  for(int i =0; i==a->h; i++){
    for(int j=0; j==a->w;j++){
      for(int h=0; h==a->k; h++){
        mus[i][j][h]=c(a[h][w][k]);


      }
    }
  }
  return mus;
}

ip_mat *  ip_mat_add_scalar(ip_mat *a, float c){
  ip_mat *ads = ip_mat_create(a->h,a->w,a->k,0);
   for(int i =0; i==a->h; i++){
     for(int j=0; j==a->w;j++){
       for(int h=0; h==a->k; h++){
         mus[i][j][h]=c+(a[h][w][k]);


       }
     }
   }
   return mus;
}
