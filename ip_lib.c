/*
 Created by Sebastiano Vascon on 23/03/20.
*/
/* id : 75; Alessio Campanelli 878170, Lorenzo Vigoni 880299, Mario Gottardo 879088, Alberto Baesso 880111*/
#include <stdio.h>
#include "ip_lib.h"
#include "bmp.h"


ip_mat * ip_mat_create(unsigned int h, unsigned int w, unsigned int k, float v){
  unsigned int i,j,l;
  ip_mat *nuova;
  stats * st;
  float ***p3;
  if (h<=0 || w<=0 || k<= 0) {
    printf("Dimensioni non valide\n");
    exit(1);
  }

  nuova = malloc(sizeof(ip_mat));
  nuova->w = w;
  nuova->h = h;
  nuova->k = k;

  st = malloc(k * sizeof(stats));
  for ( i = 0; i < k; i++){
    st[i].min =  v;
    st[i].max =  v;
    st[i].mean = v;
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
  unsigned int i,j;
  if(t){
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
}


float get_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k){
    if(i<a->h && j<a->w &&k<a->k){  /* j>=0 and k>=0 and i>=0 is non sense*/
        return a->data[i][j][k];
    }else{
        printf("Errore get_val!!!\n");
        exit(1);
    }
}


void set_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k, float v){
    if(i<a->h && j<a->w &&k<a->k){
        a->data[i][j][k]=v;
    }else{
        printf("Errore set_val!!!\n");
        exit(1);
    }
}


void compute_stats(ip_mat * t){
  unsigned int i, j, l;
  float min,max,mean;

  if( !t ){
    printf("Matrice non valida\n");
    exit(1);
  }

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

    t->stat[l].min = min;
    t->stat[l].max = max;
    t->stat[l].mean = mean;
  }
}


void ip_mat_init_random(ip_mat * t, float mean, float std){
  unsigned int i, j, l;
  if( !t ){
    printf("Matrice non valida\n");
    exit(1);
  }

  for (i = 0; i < t->h; i++){
    for (j = 0; j < t->w; j++){
      for (l = 0; l < t->k; l++){
        t->data[i][j][l] = get_normal_random(mean, std);
      }
    }
  }
  compute_stats(t);
}


ip_mat * ip_mat_copy(ip_mat * in){
    unsigned int i,j,l;
    ip_mat * copia;
    if( !in ){
      printf("Matrice non valida\n");
      exit(1);
    }

    copia = ip_mat_create(in->h, in->w, in->k, 0.0);
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
  unsigned int i, j, l;
  unsigned int row = row_end-row_start;
  unsigned int col = col_end-col_start;
  ip_mat * nuova;

  if(row <= 0 || col <= 0 || !t ||row_end > t->h || col_end > t->w ){
    printf("Matrice non valida o valori errati\n");
    exit(1);
  }

  nuova =  ip_mat_create(row, col, t->k, 0.0);

  for(i = 0; i < row; i++){
    for(j = 0; j < col; j++){
      for(l = 0; l < t->k; l++){
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
  unsigned int di=0, dj=0, dl=0;

  if(!a || !b){
    printf("Matrice non valida\n");
    exit(1);
  }

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
    printf("Almeno due dimensioni su tre non coincidono\n");
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
  unsigned int i, j, l;
  ip_mat * sum;

  if (!a || !b || a->w != b->w || a->h != b->h || a->k != b->k){
    printf("Matrice nulla o dimensioni non coincidenti\n");
    exit(1);
  }

  sum = ip_mat_copy(a);
  for (i = 0; i < a->h; i++) {
    for(j = 0; j < a->w; j++){
      for(l = 0; l < a->k; l++){
        sum->data[i][j][l] += b->data[i][j][l];
      }
    }
  }
  compute_stats(sum);
  return sum;
}


ip_mat * ip_mat_sub(ip_mat * a, ip_mat * b){
    unsigned int i, j, l;
    ip_mat * sub;

    if (!a || !b || a->w != b->w || a->h != b->h || a->k != b->k){
      printf("Matrice nulla o dimensioni non coincidenti\n");
      exit(1);
    }

    sub = ip_mat_copy(a);
    for (i = 0; i < a->h; i++) {
        for(j = 0; j < a->w; j++){
            for(l = 0; l < a->k; l++){
                sub->data[i][j][l] -= b->data[i][j][l];
            }
        }
    }
    compute_stats(sub);
    return sub;
}


ip_mat * ip_mat_mul_scalar(ip_mat *a, float c){
  unsigned int i, j, l;
  ip_mat * mus;

  if (!a){
    printf("Matrice non valida\n");
    exit(1);
  }

  mus = ip_mat_copy(a);
  for( i=0; i<a->h; i++){
    for( j=0; j<a->w;j++){
      for( l=0; l<a->k; l++){
        mus->data[i][j][l] *= c;
      }
    }
  }

  compute_stats(mus);
  return mus;
}


ip_mat *  ip_mat_add_scalar(ip_mat *a, float c){
  unsigned int i, j, l;
  ip_mat * ads;

  if (!a){
    printf("Matrice non valida\n");
    exit(1);
  }

  ads = ip_mat_copy(a);
  for( i=0; i<a->h; i++){
    for( j=0; j<a->w;j++){
      for( l=0; l<a->k; l++){
        ads->data[i][j][l] += c;
      }
    }
  }

  compute_stats(ads);
  return ads;
}

/* Moltiplica un ip_mat "a" per un ip_mat "b". Il prodotto tra i due tensori avviene cella per cella. 
 *
 * c[i][j][k] = a[i][j][k] * b[i][j][k]
 *
 * Se "a" e "b" hanno dimensioni diverse allora l'operazione non e' possibile ( uscite con exit(1)).
 *
 * Il risultato viene salvato e restituito in output all'interno di una nuova ip_mat.
 * 
 */
ip_mat * ip_mat_mul(ip_mat *a, ip_mat * b){
	unsigned int i, j, l;
	ip_mat * mul = ip_mat_create(a->h,a->w,a->k,0.0);
	if (a->h != b->h || a->w != b->w || a->k != b->k)
	{
		printf("L'operazione non e' possibile");
		exit(1);
	}
	for (int i = 0; i < a->h; i++)
	{
		for (int j = 0; j < a->w; j++)
		{
			for (int l = 0; l < a->k; l++)
			{
				mul->data[i][j][l] = a->data[i][j][l] * b->data[i][j][l];
			}
		}
	}
	return mul;
}

ip_mat * ip_mat_mean(ip_mat * a, ip_mat * b){
  unsigned int i, j, l;
  ip_mat * mean;

  if (!a || !b || a->w != b->w || a->h != b->h || a->k != b->k){
    printf("Matrice nulla o dimensioni non coincidenti\n");
    exit(1);
  }

  mean = ip_mat_copy( a );
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

ip_mat * ip_mat_to_gray_scale(ip_mat * in){
    unsigned int i, j, l;
    ip_mat *bw;

    if( !in ){
      printf("Matrice non valida\n");
      exit(1);
    }

    bw = ip_mat_create(in->h, in->w, in->k, 1.0);
    for(i = 0; i < in->h; i++){
        for(j = 0; j < in->w; j++){
          float mean = 0;
            for(l = 0; l < in->k; l++){
              mean += in->data[i][j][l] / in->k;
            }
            for(l = 0; l < in->k; l++){
              bw->data[i][j][l] = mean;
            }
        }
    }
    compute_stats(bw);
    return bw;
}


ip_mat * ip_mat_warhol(ip_mat * in){
	unsigned int i, j;
	ip_mat * switchRG, * switchGB, * switchRB;/*lavoro sul k della matrice*/
	ip_mat * out, *half_left, *half_right;

	switchRG = ip_mat_create(in->h,in->w,in->k,0);
	switchGB = ip_mat_create(in->h,in->w,in->k,0);
	switchRB = ip_mat_create(in->h,in->w,in->k,0);

	for (i = 0; i < in->h; i++) {
    	for(j = 0; j < in->w; j++){
        	switchRG->data[i][j][0] = in->data[i][j][1]; /*top dx*/
            switchRG->data[i][j][1] = in->data[i][j][0];
            switchRG->data[i][j][2] = in->data[i][j][2];

            switchGB->data[i][j][1] = in->data[i][j][2]; /*bottom sx*/
            switchGB->data[i][j][2] = in->data[i][j][1];
			switchGB->data[i][j][0] = in->data[i][j][0];

            switchRG->data[i][j][0] = in->data[i][j][2]; /*bottom dx*/
            switchRG->data[i][j][2] = in->data[i][j][0];
            switchRG->data[i][j][1] = in->data[i][j][1];
      	}
  	}
  	half_left = ip_mat_concat(in,switchGB,0);
  	half_right = ip_mat_concat(switchRG,switchRB,0);
  	out = ip_mat_concat(half_left,half_right,1);
  	ip_mat_free(switchRG):
  	ip_mat_free(switchGB);
  	ip_mat_free(switchRB):
  	ip_mat_free(half_left);
  	ip_mat_free(half_right):

  	return out;

}

ip_mat * ip_mat_blend(ip_mat * a, ip_mat * b, float alpha) {
  ip_mat * blended;
  unsigned int i,j,l;

  if(!a || !b || a->h != b->h || a->w != b->w || a->k != b->k || alpha < 0 || alpha > 1){
    printf("Immagini di dimensione diversa o alpha di dimensione errata\n");
    exit(1);
  }

  blended = ip_mat_create(a->h,a->w,a->k,0);
  for (i = 0; i < blended->h; i++) {
    for(j = 0; j < blended->w; j++){
      for(l = 0; l < blended->k; l++){
        blended->data[i][j][l] =  (alpha *a->data[i][j][l]) + ((1 - alpha) * b->data[i][j][l]);
      }
    }
  }
  return blended;
}


ip_mat * ip_mat_brighten(ip_mat * a, float bright){
  ip_mat * lux;

  if(!a){
    printf("Matrice non valida\n");
    exit(1);
  }

  lux = ip_mat_add_scalar(a,  bright);
  compute_stats (lux);
  return lux;
}

ip_mat * ip_mat_corrupt( ip_mat * a, float amount ){
    ip_mat * b;

    if (!a || amount < 0 || amount > 255){
      printf("Amount non compreso tra 0 e 255 o matrice nulla\n");
      exit(1);
    }

    b = ip_mat_copy(a);
    ip_mat_init_random(b, 0, amount/2);
    b = ip_mat_sum(a, b);
    compute_stats(b);
    return b;
}

ip_mat * ip_mat_convolve(ip_mat * a, ip_mat * f){
  int deltaW, deltaH;
  ip_mat * filtered, *out = ip_mat_create(a->h, a->w, a->k, 0);
  ip_mat * filter;
  unsigned int i, j, I, J, L;

  deltaW = (f->w -1) / 2;
  deltaH = (f->h -1) / 2;

  if( a->k % f->k == 0 && a && f){
    filter = ip_mat_copy(f);
    while(a->k > filter->k){
      filter = ip_mat_concat(filter, f, 2);
    }
  }
  else{
    printf("Dimensioni filtro/foto non compatibili o matrie nulla\n");
    exit(1);
  }

  filtered = ip_mat_padding(a, deltaH, deltaW );

  for( I = deltaH; I < filtered->h - deltaH; I++ ){
  	for( J = deltaW; J < filtered->w - deltaW; J++ ){
  	  for( L = 0; L < filtered->k; L++ ){

  	  	for( i = 0; i < filter-> h; i++ ){
  	  	  for( j= 0; j < filter-> w; j++ ){
  	  	  	out->data[I-deltaH][J-deltaW][L] += filter->data[i][j][L] *filtered->data[I-deltaH+i][J-deltaW+j][L];
  	  	  }
  	  	}

  	  }
  	}
  }
  ip_mat_free(filter);
  ip_mat_free(filtered);
  compute_stats(out);
  return out;
}


ip_mat * ip_mat_padding(ip_mat * a, unsigned int pad_h, unsigned int pad_w){
  ip_mat *padd;
  unsigned int i,j,l;

  if (!a){
    printf("Matrice non valida\n");
    exit(1);
  }

  padd = ip_mat_create (a->h+(pad_h*2), a->w+(pad_w*2), a->k, 0);

  for (i = 0; i < a->h; i++) {
      for(j = 0; j < a->w; j++){
          for(l = 0; l < padd->k; l++){
              padd->data[i+pad_h][j+pad_w][l] = a->data[i][j][l];
          }
      }
  }
  compute_stats(padd);
  return padd;
}


ip_mat * create_sharpen_filter(){
  ip_mat *sharpen = ip_mat_create(3,3,1,0.0);
  set_val(sharpen,0,1,0,-1);
  set_val(sharpen,1,0,0,-1);
  set_val(sharpen,1,2,0,-1);
  set_val(sharpen,2,1,0,-1);
  set_val(sharpen,1,1,0,5);
  compute_stats(sharpen);
  return sharpen;
}

ip_mat * create_edge_filter(){
  ip_mat *edge = ip_mat_create(3,3,1,-1.0);
  set_val(edge,1,1,0,8.0);
  compute_stats(edge);
  return edge;
}

ip_mat * create_sobel_horizontal(){
  ip_mat *sobelHzl = ip_mat_create(3,3,1,0.0);
  set_val(sobelHzl,0,0,0,-1.0);
  set_val(sobelHzl,0,2,0,-1.0);
  set_val(sobelHzl,0,1,0,2.0);
  set_val(sobelHzl,2,1,0,2.0);
  set_val(sobelHzl,2,0,0,1.0);
  set_val(sobelHzl,2,2,0,1.0);
  compute_stats(sobelHzl);
  return sobelHzl;	
}

ip_mat * create_sobel_vertical(){
  ip_mat *sobelVtl = ip_mat_create(3,3,1,0.0);
  set_val(sobelVtl,0,0,0,-1.0);
  set_val(sobelVtl,2,0,0,-1.0);
  set_val(sobelVtl,1,0,0,2.0);
  set_val(sobelVtl,1,2,0,2.0);
  set_val(sobelVtl,0,2,0,1.0);
  set_val(sobelVtl,2,2,0,1.0);
  compute_stats(sobelVtl);
  return sobelVtl;
}

ip_mat * create_emboss_filter(){
  ip_mat *emboss = ip_mat_create(3,3,1,1.0);
  set_val(emboss,0,0,0,-2.0);
  set_val(emboss,0,1,0,-1.0);
  set_val(emboss,0,2,0,0.0);
  set_val(emboss,1,0,0,-1.0);
  set_val(emboss,2,0,0,0.0);
  set_val(emboss,2,2,0,2.0);
  compute_stats(emboss);
  return emboss;
}

ip_mat * create_average_filter(unsigned int w, unsigned int h, unsigned int k){
  float c=1.0/(w*h);

  ip_mat * avg;

  if (h%2 == 0 || w%2 == 0 || k == 0){
    printf("Il filtro deve avere dimensioni dispari e k > 0\n");
    exit(1);
  }

  avg = ip_mat_create(w,h,k,c);
  compute_stats(avg);
  return avg;
}

ip_mat * create_gaussian_filter(unsigned int h, unsigned int w, unsigned int k, float sigma){
    ip_mat * gaussian;
    unsigned int i, j, l;
    int x, y, cx, cy;
    float sum = 0.0;

    if (h%2 == 0 || w%2 == 0 || k == 0){
        printf("Il filtro deve avere dimensioni dispari e k > 0\n");
        exit(1);
    }

    gaussian = ip_mat_create(h, w, k, 0.0);
    cx = (w-1) / 2;
    cy = (h-1) / 2;
    for(i = 0; i < h; i++){
        for(j = 0; j < w; j++){
            float value;
            y = i - cy;
            x = j - cx;
            value = (1.0/(2.0*PI*sigma*sigma))*exp(-(x*x+y*y)/(2.0*sigma*sigma));
            sum += value;
            for(l = 0; l < k; l++){
              set_val(gaussian, i, j, l, value);
            }
        }
    }
    gaussian = ip_mat_mul_scalar(gaussian, 1/sum);

    compute_stats(gaussian);
    return gaussian;
}

void rescale(ip_mat * t, float new_max){
  unsigned int i, j, l;

  if(!t){
    printf("Matrice non valida\n");
    exit(1);
  }

  for(i=0; i<t->h; i++){
    for(j=0; j<t->w; j++){
      for(l=0; l<t->k; l++){
        t->data[i][j][l] = (t->data[i][j][l] - t->stat[l].min)/(t->stat[l].max - t->stat[l].min) * new_max;
      }
    }
  }
  compute_stats(t);
}

void clamp(ip_mat * t, float low, float high){
  unsigned int i, j, l;

  if(!t){
    printf("Matrice non valida\n");
    exit(1);
  }

  for(i=0; i<t->h; i++){
    for(j=0; j<t->w; j++){
      for(l=0; l<t->k; l++){
        if (t->data[i][j][l] > high)
          t->data[i][j][l] = high;
        if (t->data[i][j][l] < low)
          t->data[i][j][l] = low;
      }
    }
  }
  compute_stats(t);
}

float get_normal_random(float media, float std){

    float y1 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    float y2 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    float num = cos(2*PI*y2)*sqrt(-2.*log(y1));

    return media + num*std;
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
