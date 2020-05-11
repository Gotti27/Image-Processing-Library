/*
 Created by Sebastiano Vascon on 23/03/20.
*/

#include <stdio.h>
#include "ip_lib.h"
#include "bmp.h"


ip_mat * ip_mat_create(unsigned int h, unsigned int w, unsigned int k, float v){
  unsigned int i,j,l;
  ip_mat *nuova;
  stats * st;
  float ***p3;

  nuova = malloc(sizeof(ip_mat));
  nuova->w = w;
  nuova->h = h;
  nuova->k = k;

  st = malloc(k * sizeof(stats));
  for ( i = 0; i < k; i++){
    (&st[i])->min =  v;
    (&st[i])->max =  v;
    (&st[i])->mean = v;
  }
  nuova->stat = st;

  p3 = malloc (h * sizeof(float**));
  for ( i = 0; i < h; i++){
    p3[i] =  malloc (w * sizeof(float*));
    for ( j = 0; j < w; j++){
      p3[i][j] = malloc (k * sizeof(float));
      for ( l = 0; l < k; l++){
        p3[i][j][l] = v;
      }
    }
  }
  nuova->data = p3;

  return nuova;
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


void compute_stats(ip_mat * t){
  unsigned int i, j, l;
  float min,max,mean;

  for (l = 0; l < t->k; l++){
    min =  t->data[0][0][l];
    max =  t->data[0][0][l];
    mean = 0.0;

    for ( i = 0; i < t->h; i++){
      for ( j = 0; j < t->w; j++){
        if(t->data[i][j][l] > max ){
          max = t->data[i][j][l];
        }
        if(t->data[i][j][l] < min ){
          min = t->data[i][j][l];
        }
        mean += t->data[i][j][l] / (t->h * t->w);
      }
    }

    (&t->stat[l])->min = min;
    (&t->stat[l])->max = max;
    (&t->stat[l])->mean = mean;
  }
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


ip_mat * ip_mat_subset(ip_mat * t, unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end){
  int i, j, l;
  int rig=row_end-row_start;
  int col= col_end-col_start;
  ip_mat * nuova =  ip_mat_create(rig, col, t->k, 0.0);

  for(i=0; i<rig; i++){
    for(j=0; j<col; j++){
      for(l=0; l<t->k; l++){
        nuova->data[i][j][l] = t->data[i+row_start][j+col_start][l];
      }
    }
  }

  compute_stats(nuova);
  return nuova;
}


ip_mat * ip_mat_concat(ip_mat * a, ip_mat * b, int dimensione){
  ip_mat * out;
  unsigned int i,j,l;
  int di=0, dj=0, dl=0;

  if (dimensione == 0 && a->w == b-> w && a->k == b->k){
    out = ip_mat_create((a->h + b->h), a->w,a->k,  0);
    di = a->h;
  }
  else if(dimensione == 1 && a->h == b->h && a->k == b->k){
    out = ip_mat_create(a->h, (a->w + b->w), a->k, 0);
    dj = a->w;
  }
  else if(dimensione == 2 && a->h == b->h && a->w == b->w){
    out = ip_mat_create(a->h, a->w, (a->k + b->k),0);
    dl = a->k;
  }
  else{
    exit(1);
  }

  for (i = 0; i < out->h; i++){
    for (j = 0; j < out->w; j++){
      for (l = 0; l < out->k; l++){
        if (i<di || j<dj || l<dl){
          out->data[i][j][l] = a->data[i][j][l];
        } else {
          out->data[i][j][l] = b->data[i-di][j-dj][l-dl];
        }
      }
    }
  }

  compute_stats(out);
  return out;
}

ip_mat * ip_mat_sum(ip_mat * a, ip_mat * b){
  if (a->w == b->w && a->h == b->h && a->k == b->k) {
      int i, j, l;
      ip_mat * sum = ip_mat_copy(a);

      for (i = 0; i < a->h; i++) {
          for(j = 0; j < a->w; j++){
              for(l = 0; l < a->k; l++){
                  sum->data[i][j][l] = a->data[i][j][l] + b->data[i][j][l];
              }
          }
      }
      compute_stats(sum);
      return sum;
    }
  else{
    exit(1);
  }
}


ip_mat * ip_mat_sub(ip_mat * a, ip_mat * b){
    if (a->w == b->w && a->h == b->h && a->k == b->k) {
        int i, j, l;
        ip_mat * sub = ip_mat_copy(a);

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
  int i, j, l;
  ip_mat * mus = ip_mat_copy(a);

  for( i=0; i<a->h; i++){
    for( j=0; j<a->w;j++){
      for( l=0; l<a->k; l++){
        mus->data[i][j][l]=c*(a->data[i][j][l]);
      }
    }
  }

  compute_stats(mus);
  return mus;
}


ip_mat *  ip_mat_add_scalar(ip_mat *a, float c){
  int i, j, l;
  ip_mat * ads = ip_mat_copy(a);

  for( i=0; i<a->h; i++){
    for( j=0; j<a->w;j++){
      for( l=0; l<a->k; l++){
        ads->data[i][j][l]=c+(a->data[i][j][l]);
      }
    }
  }

  compute_stats(ads);
  return ads;
}

ip_mat * ip_mat_mean(ip_mat * a, ip_mat * b){
    if(a->w == b->w && a->h == b->h && a->k == b->k){
        int i, j, l;
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

ip_mat * ip_mat_to_gray_scale(ip_mat * in){
    int i, j, l;
    ip_mat *bw = ip_mat_create(in->h, in->w, in->k, 1.0);
    for(i = 0; i < in->h; i++){
        for(j = 0; j < in->w; j++){
            int mean = (in->data[i][j][0] + in->data[i][j][1] + in->data[i][j][2]) / 3;
            for(l = 0; l < in->k; l++){
                bw->data[i][j][l] = mean;
            }
        }
    }
    compute_stats(bw);
    return bw;
}

/*
ip_mat * ip_mat_blend(ip_mat * a, ip_mat * b, float alpha) {
  ip mat * blended;
  if(a->h == b->h && a->w == b->w && a->k == b->k && alpha >=0 && alpha <=1)
    {
      unsigned int i,j,l;
      blended = ip_mat_create(a->h,a->w,a->k,0);
      for (i = 0; i < blended->h; i++) {
            for(j = 0; j < blended->w; j++){
                for(l = 0; l < blended->k; l++){
                    blended->data[i][j][l] =  (alpha *a->data[i][j][l]) + ((1 - alpha) * b->data[i][j][l]) ;
                }
            }
        }
    }
    else
    {
      printf("immagini di dimensione diversa o alpha di dimensione errata\n");
      exit(1);

    }
}
*/
ip_mat * ip_mat_blend(ip_mat * a, ip_mat * b, float alpha) { /* V2.0 ---- da due immagini ne da un'altra creata dalle dimensioni minori*/
  ip_mat * blended;
  if(alpha >=0 && alpha <=1)
    {
      unsigned int i,j,l,h,w,k;
      h = a->h;
      w = a->w;
      k = a->k;
      blended = ip_mat_create(a->h,a->w,a->k,0);
      for (i = 0; i < blended->h; i++) {
            for(j = 0; j < blended->w; j++){
                for(l = 0; l < blended->k; l++){
                  if (i < h && j< w && k < k) {
                    blended->data[i][j][l] =  (alpha *a->data[i][j][l]) + ((1 - alpha) * b->data[i][j][l]) ;
                  }
                  else {
                    blended->data[i][j][l] = a->data[i][j][l];
                  }
                }
            }
        }
    }
    else
    {
      printf(" alpha di dimensione errata\n");
      exit(1);

    }
    compute_stats(blended);
    return blended;
}

ip_mat * ip_mat_brighten(ip_mat * a, float bright){
 ip_mat * lux = ip_mat_add_scalar( a,  bright);
 clamp(lux,0.0,255.0);
 compute_stats (lux);
 return lux;
}

ip_mat * ip_mat_corrupt( ip_mat * a, float amount ){
 if (amount < 0 || amount > 255){
   exit(1);
 }
  ip_mat * b = ip_mat_copy(a);
  int mean = ((&a->stat[0])->mean + (&a->stat[1])->mean + (&a->stat[2])->mean ) /3;

  ip_mat_init_random(b, mean, amount);
  b = ip_mat_sum(a, b);
  compute_stats(b)

  return b;
}

ip_mat * ip_mat_convolve(ip_mat * a, ip_mat * f){
  int deltaW, deltaH;
  ip_mat * filtered;
  int i, j, l, I, J, L;

  deltaW = (f->w -1) / 2;
  deltaH = (f->h -1) / 2;

  filtered = ip_mat_create( a->h - deltaH * 2, a->w - deltaW * 2, a->k, 0);

  for( I = deltaH; I < filtered->h - deltaH; I++ ){
  	for( J = deltaW; J < filtered->w - deltaW; J++ ){
  	  for( L = 0; L < filtered->k; L++ ){

  	  	for( i = 0; i < f->h; i++ ){
  	  	  for( j= 0; j < f-> w; j++ ){
  	  	  	filtered->data[I][J][L] += f->data[i][j][0] * a->data[I-deltaH+i][J-deltaW+j][L];
  	  	  }
  	  	}

  	  }
  	}
  }
  return filtered;
}

void rescale(ip_mat * t, float new_max){
  int i, j, l;

  for(i=0; i<t->h; i++){
    for(j=0; j<t->w; j++){
      for(l=0; l<t->k; l++){
        t->data[i][j][l] = (t->data[i][j][l] - (&t->stat[l])->min)/((&t->stat[l])->max - (&t->stat[l])->min) * new_max;
      }
    }
  }
}

void clamp(ip_mat * t, float low, float high){
  int i, j, l;

  for(i=0; i<t->h; i++){
    for(j=0; j<t->w; j++){
      for(l=0; l<t->k; l++){
        if (t->data[i][j][l] > high)
          t->data[i][j][l] = 255.0;
        if (t->data[i][j][l] > high)
          t->data[i][j][l] = 0.0;
      }
    }
  }
}

float get_normal_random(){
    float y1 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    float y2 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    return cos(2*PI*y2)*sqrt(-2.*log(y1));

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
