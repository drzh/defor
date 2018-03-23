/*
 * Usage: prog -i FILE_chr_pos_freq(default: 'stdin') -w INT_window(1000000) -f FLOAT_min_mean_cutoff(0.05) -F FLOAT_max_mean_cutoff(0.95) -I INT_max_iteration | <STDOUT>
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <gsl/gsl_randist.h>

typedef struct Param_Struct {
  char file1[1000];
  char file2[1000];
  long w;
  double f;
  double F;
  int I;
} Param;

typedef struct Info_Struct {
  long id;
  char * chr;
  long pos;
  double frq[2];
  struct Info_Struct * next;
} Info;
  
FILE * xfopen (char * file_name, char * p);
void * xfread (char * file_name);
void xfwrite (void * p, int size, char * file_name);
void * xmalloc(size_t size);
void * xrealloc(void *ptr, size_t size);
void xfree(void *ptr);
void usage(char * prog);
void die(const char *fmt, ...);
Param * getoption(int argc, char *argv[]);
double mean(double * x, double * prob, int n);
double sd(double * x, double * prob, int n);
void nclust(Info * start, Info * end, int paramI, int idxn);
void readinfo(char * line, Info * info);

int main (int argc, char * argv[]) {
  // Read data
  Param * param = getoption(argc, argv);
  Info * infostart = NULL;
  Info * infoend = NULL;
  int infosize = sizeof(Info);
  char * line1 = NULL;
  char * line2 = NULL;
  long n1 = 0;
  long n2 = 0;
  Info * info1 = xmalloc(infosize);
  Info * info2 = xmalloc(infosize);
  info1 -> chr = xmalloc(1000);
  info2 -> chr = xmalloc(1000);
  FILE * fp1 = NULL;
  FILE * fp2 = NULL;
  int nidx = 0;
  if (param -> file1[0]) {
    nidx = 1;
    fp1 = xfopen(param -> file1, "r");
  }
  else {
    usage(argv[0]);
  }
  if (param -> file2[0]) {
    nidx = 2;
    fp2 = xfopen(param -> file2, "r");
  }
  while (getline(&line1, &n1, fp1) != -1) {
    readinfo(line1, info1);
    while (nidx == 2 && (! line2 || strverscmp(info1 -> chr, info2 -> chr) > 0 || strverscmp(info1 -> chr, info2 -> chr) == 0 && info1 -> pos > info2 -> pos) && getline(&line2, &n2, fp2) != -1) {
      readinfo(line2, info2);
    }
    if ((nidx == 1 &&
    	 info1 -> frq[0] >= param -> f && info1 -> frq[0] <= param -> F
    	 ) ||
    	(nidx == 2 &&
    	 strverscmp(info1 -> chr, info2 -> chr) == 0 &&
    	 info1 -> pos == info2 -> pos &&
    	 (info1 -> frq[0] >= param -> f && info1 -> frq[0] <= param -> F ||
    	  info2 -> frq[0] >= param -> f && info2 -> frq[0] <= param -> F
    	  )
    	 )
    	) {
      Info * info = xmalloc(infosize);
      memset(info, 0, infosize);
      info -> chr = xmalloc((strlen(info1 -> chr) + 1) * sizeof(char));
      strcpy(info -> chr, info1 -> chr);
      info -> pos = info1 -> pos;
      info -> frq[0] = info1 -> frq[0];
      info -> frq[1] = info2 -> frq[0];
      if (infoend) {
    	info -> id = infoend -> id + 1;
    	infoend -> next = info;
    	infoend = info;
      }
      else {
    	info -> id = 1;
    	infostart = info;
    	infoend = info;
      }
    }
  }
  fclose(fp1);
  fclose(fp2);
  
  // Process data
  Info * start = infostart;
  Info * end = infostart;
  while (end) {
    while (end -> next &&
  	   strverscmp(start -> chr, end -> next -> chr) == 0 &&
  	   start -> pos + param -> w > end -> next -> pos
  	   ) {
      end = end -> next;
    }
    nclust(start, end, param -> I, nidx);
    end = end -> next;
    if (end) {
      if (strverscmp(start -> chr, end -> chr) == 0) {
  	while (start -> pos + param -> w < end -> pos) {
  	  start = start -> next;
  	}
      }
      else {
  	start = end;
      }
    }
  }
  // Free memory
  Info * info = infostart;
  while (info) {
    Info * p = info;
    info = info -> next;
    xfree(p -> chr);
    xfree(p);
  }
  xfree(info1 -> chr);
  xfree(info1);
  xfree(info2 -> chr);
  xfree(info2);
}

void readinfo(char * line, Info * info) {
  char delimiters[] = "\t\n";
  char * ele = NULL;
  char * running = line;
  int i = 0;
  while ((ele = strsep(&running, delimiters)) && ++i < 4) {
    if (i == 1) {
      strcpy(info -> chr, ele);
    }
    else if (i == 2) {
      char * tail;
      info -> pos = strtol(ele, &tail, 0);
    }
    else if (i == 3) {
      char * tail;
      info -> frq[0] = strtod(ele, &tail);
    }
  }
  return;
}

void nclust(Info * start, Info * end, int paramI, int nidx) {
  long n = end -> id - start -> id + 1;
  int k = 4;
  double * x = xmalloc(n * sizeof(double));
  double * xlhsum = xmalloc(n * sizeof(double));
  double * xmean = xmalloc(k * sizeof(double));
  double * xmeanpre = xmalloc(k * sizeof(double));
  double * xsd = xmalloc(k * sizeof(double));
  double * xlh = xmalloc(n * k * sizeof(double));
  double * xprob = xmalloc(n * k * sizeof(double));
  int * xgrp = xmalloc(n * sizeof(int));
  long * xgrpcnt = xmalloc(k * sizeof(long));
  /* double * xgrpprop = xmalloc(k * sizeof(double)); */
  int i, j, idxgrp;
  printf("%s\t%d\t%d", start -> chr, start -> pos, end -> pos);
  for (idxgrp = 0; idxgrp < nidx; idxgrp++) {
    int flag = 1, ite = 0;
    Info * info = start;
    for (i = 0; i < n; i++) {
      x[i] = info -> frq[idxgrp];
      info = info -> next;
    }
    // Initialize parameters
    for (i = 0; i < k; i++) {
      xmean[i] = (double) i / (k - 1);
      xsd[i] = 0.05;
    }
    for (i = 0; i < n * k; i++) {
      xprob[i] = 0.25;
    }
    while (flag && ite++ < paramI) {
      memset(xlhsum, 0, n * sizeof(double));
      for (i = 0; i < k; i++) {
	for (j = 0; j < n; j++) {
	  /* xlh[i * n + j] = xsd[i] ? gsl_ran_gaussian_pdf(x[j] - xmean[i], xsd[i]) : 0; */
	  // Calc likelihood of point within 2 * SD of mean
	  xlh[i * n + j] = (xmean[i] - 2 * xsd[i] < x[j] && xmean[i] + 2 * xsd[i] > x[j]) ? gsl_ran_gaussian_pdf(x[j] - xmean[i], xsd[i]) : 0;
	  xlhsum[j] += xlh[i * n + j];
	}
      }
      for (i = 0; i < k; i++) {
	for (j = 0; j < n; j++) {
	  xprob[i * n + j] = xlhsum[j] > 0 ? xlh[i * n + j] / xlhsum[j] : 0;
	}
      }
      memcpy(xmeanpre, xmean, k * sizeof(double));
      flag = 0;
      for (i = 0; i < k; i++) {
	xmean[i] = mean(x, xprob + i * n, n);
	xsd[i] = sd(x, xprob + i * n, n);
	if (xmean[i] - xmeanpre[i] > 0.001) {
	  flag = 1;
	}
      }
    }
    memset(xgrp, 0, n * sizeof(int));
    for (i = 1; i < k; i++) {
      for (j = 0; j < n; j++) {
	if (xprob[xgrp[j] * n + j] < xprob[i * n + j]) {
	  xgrp[j] = i;
	}
      }
    }
    memset(xgrpcnt, 0, k * sizeof(long));
    for (j = 0; j < n; j++) {
      xgrpcnt[xgrp[j]]++;
    }
    /* memset(xgrpprop, 0, k * sizeof(double)); */
    /* for (i = 0; i < k; i++) { */
    /*   xgrpprop[i] = (double) xgrpcnt[i] / n; */
    /* } */
  // Output result
  /* printf("%s\t%d\t%d\t", start -> chr, start -> pos, end -> pos); */
    printf("\t");
    for (i = 0; i < k; i++) {
      if (i) {
	printf(",");
      }
      printf("%.3f", xmean[i]);
    }
    printf("\t");
    for (i = 0; i < k; i++) {
      if (i) {
	printf(",");
      }
      printf("%.3f", xsd[i]);
    }
    printf("\t");
    for (i = 0; i < k; i++) {
      if (i) {
	printf(",");
      }
      printf("%d", xgrpcnt[i]);
    }
    /* printf("\t"); */
    /* for (i = 0; i < k; i++) { */
    /*   if (i) { */
    /* 	printf(","); */
    /*   } */
    /*   printf("%.3f", xgrpprop[i]); */
    /* } */
  }
  printf("\n");
  xfree(x);
  xfree(xlhsum);
  xfree(xmean);
  xfree(xmeanpre);
  xfree(xsd);
  xfree(xlh);
  xfree(xprob);
  xfree(xgrp);
  xfree(xgrpcnt);
  /* xfree(xgrpprop); */
  return;
}

double mean(double * x, double * prob, int n) {
  double sum = 0;
  double nprob = 0;
  int i;
  for (i = 0; i < n; i++) {
    sum += prob[i] * x[i];
    nprob += prob[i];
  }
  return nprob > 0 ? sum / nprob : 0;
}

double sd(double * x, double * prob, int n) 
{
  double a = 0;
  double b = 0;
  double xmean = mean(x, prob, n);
  int i;
  for (i = 0; i < n; i++) {
    a += prob[i] * (x[i] - xmean) * (x[i] - xmean);
    b += prob[i];
  }
  return b > 0 ? sqrt(a / b) : 0;
}

Param * getoption(int argc, char *argv[]){
  Param * param = xmalloc(sizeof(Param));
  param -> w = 10000000;
  param -> f = 0;
  param -> F = 1;
  param -> I = 10;
  memset(param -> file1, 0, 1000);
  memset(param -> file2, 0, 1000);
  char c;
  extern char *optarg;
  extern int optind, optopt;
  while ((c = getopt(argc, argv, "i:w:n:d:f:F:p:c:t:I:h")) != (char)-1) {
    switch(c) {
      /* case 'i': */
      /*   if (sscanf(optarg, "%s", param -> filei) != 1) */
      /*     usage(argv[0]); */
      /*   break; */
      case 'w':
        if (sscanf(optarg, "%ld", &param -> w) != 1)
          usage(argv[0]);
        break;
      case 'f':
        if (sscanf(optarg, "%lf", &param -> f) != 1)
          usage(argv[0]);
        break;
      case 'F':
        if (sscanf(optarg, "%lf", &param -> F) != 1)
          usage(argv[0]);
        break;
      case 'I':
        if (sscanf(optarg, "%d", &param -> I) != 1)
          usage(argv[0]);
        break;
    }
  }
  if (optind < argc) {
    strcpy(param -> file1, argv[optind++]);
  }
  else {
    usage(argv[0]);
  }
  if (optind < argc) {
    strcpy(param -> file2, argv[optind++]);
  }
  return param;
}

void usage(char * prog)
{
  printf("\n\
Usage: %s [options] file_normal file_tumor            \n\n\
         Options:                                    \n\
           -w : window size [1000000]                \n\
           -f : minimum allele frequency [0.05]      \n\
           -F : maximum allele frequency [0.95]      \n\
           -I : maximum iteration [10]\n             \n\
", prog);
  exit(1);
}

FILE * xfopen (char * file_name, char * p) 
{
  FILE * fp = fopen(file_name, p);
  if (! fp)
    die("Can not open file: %s\n", file_name);
  return fp;
}

void * xfread (char * file_name)
{
  FILE * fp = xfopen(file_name, "r");
  struct stat info;
  if (stat64(file_name, &info) == -1)
    die("stat64 error\n");
  void * p = xmalloc(info.st_size);
  if (fread(p, 1, info.st_size, fp)) {
    fclose(fp);
    return p;
  }
  else
    die("fread fail\n");
}

void xfwrite (void * p, int size, char * file_name) 
{
  FILE * fp = xfopen(file_name, "w");
  if (fwrite(p, 1, size, fp) != size)
    die("fwrite error\n");
  fclose(fp);
}

void * xmalloc (size_t size)
{
  void *tmp = malloc(size);
  if (tmp == NULL)
    die("Cannot malloc enough memory: %d\n", size);
  return tmp;
}

void * xrealloc(void *ptr, size_t size)
{
  void *tmp = realloc(ptr, size);
  if (tmp == NULL) 
    die("Cannot realloc enough memory!\n");
  return tmp;
}

void xfree(void *ptr)
{
  if (ptr)
    free(ptr);
  return;
}

void die(const char *fmt, ...)
{
  fprintf(stderr, "!!! Error with following message\n\t");
  va_list ap;
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
  exit(-1);
}
