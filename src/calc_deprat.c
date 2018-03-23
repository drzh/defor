/*
 * Usage: prog -w INT_window(1000000) -s INT_step(1000) -n INT_max_cluster(4) -d FLOAT_diff_cutoff(0.1) file_normal file_tumor | <STDOUT>
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
#include <gsl/gsl_sort.h>

typedef struct Param_Struct {
  char file1[1000];
  char file2[1000];
  long w;
  long s;
  int d;
  int n;
} Param;

typedef struct Info_Struct {
  char * chr;
  long pos;
  char seq;
  int dep[2];
  double ratio;
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
Info * readinfopair(FILE * fp1, FILE * fp2, Param * param);
double depratio(Info * start, long n);
double gc(Info * start, long n);

int main (int argc, char * argv[]) {
  // Read data
  Param * param = getoption(argc, argv);
  Info * start = NULL;
  Info * end = NULL;
  int infosize = sizeof(Info);
  FILE * fp1 = NULL;
  FILE * fp2 = NULL;
  if (param -> file1[0] && param -> file2[0]) {
    fp1 = xfopen(param -> file1, "r");
    fp2 = xfopen(param -> file2, "r");
  }
  else {
    usage(argv[0]);
  }
  start = readinfopair(fp1, fp2, param);
  end = start;
  long n = 1;
  long pospre = -1;
  char chrpre[1000];
  while (end) {
    Info * info = NULL;
    while ((info = readinfopair(fp1, fp2, param)) &&
	   strverscmp(start -> chr, info -> chr) == 0 &&
  	   start -> pos + param -> w > info -> pos
  	   ) {
      end -> next = info;
      end = info;
      n++;
    }
    if (n >= param -> n
	/* && */
    	/* (pospre == -1 || */
    	/*  strverscmp(chrpre, start -> chr) != 0 || */
    	/*  pospre + param -> s < start -> pos */
    	/*  ) */
    	) {
      printf("%s\t%ld\t%ld\t%.3f\t%.3f\n",
	     start -> chr,
	     start -> pos,
	     end -> pos,
	     depratio(start, n),
	     gc(start, n)
	     );
      strcpy(chrpre, start -> chr);
      pospre = start -> pos;
    }
    end -> next = info;
    end = info;
    n++;
    if (end) {
      if (strverscmp(start -> chr, end -> chr) == 0) {
  	while (start -> next &&
	       (pospre + param -> s > start -> pos
	       	||
		start -> pos + param -> w < end -> pos
		)
	       ) {
   	  Info * p = start;
  	  start = start -> next;
	  n--;
  	  xfree(p -> chr);
  	  xfree(p);
  	}
      }
      else {
  	while (start != end) {
  	  Info * p = start;
  	  start = start -> next;
	  n--;
  	  xfree(p -> chr);
  	  xfree(p);
  	}
      }
    }
  }
  while (start) {
    Info * p = start;
    start = start -> next;
    xfree(p -> chr);
    xfree(p);
  }
}

Info * readinfopair(FILE * fp1, FILE * fp2, Param * param) {
  char * line1 = NULL;
  char * line2 = NULL;
  long n1 = 0;
  long n2 = 0;
  int infosize = sizeof(Info);
  Info * info1 = xmalloc(infosize);
  Info * info2 = xmalloc(infosize);
  info1 -> chr = xmalloc(1000);
  info2 -> chr = xmalloc(1000);
  Info * info = NULL;
  while (getline(&line1, &n1, fp1) != -1) {
    readinfo(line1, info1);
    while ((! line2 || strverscmp(info1 -> chr, info2 -> chr) > 0 || strverscmp(info1 -> chr, info2 -> chr) == 0 && info1 -> pos > info2 -> pos) && getline(&line2, &n2, fp2) != -1) {
      readinfo(line2, info2);
    }
    if (strverscmp(info1 -> chr, info2 -> chr) == 0 &&
	info1 -> pos == info2 -> pos &&
	info1 -> dep[0] > param -> d &&
	info2 -> dep[0] > param -> d
	) {
      info = xmalloc(infosize);
      memset(info, 0, infosize);
      info -> chr = xmalloc(strlen(info1 -> chr) + 1);
      strcpy(info -> chr, info1 -> chr);
      info -> pos = info1 -> pos;
      info -> seq = info1 -> seq;
      info -> dep[0] = info1 -> dep[0];
      info -> dep[1] = info2 -> dep[0];
      info -> ratio = (info -> dep[0] > 0 && info -> dep[1] > 0) ? (double)info -> dep[1] / info -> dep[0] : -1;
      break;
    }
  }
  xfree(info1 -> chr);
  xfree(info2 -> chr);
  xfree(info1);
  xfree(info2);
  xfree(line1);
  xfree(line2);
  return info;
}

void readinfo(char * line, Info * info) {
  char delimiters[] = "\t\n";
  char * ele = NULL;
  char * running = line;
  int i = 0;
  while ((ele = strsep(&running, delimiters)) && ++i < 5) {
    if (i == 1) {
      strcpy(info -> chr, ele);
    }
    else if (i == 2) {
      char * tail;
      info -> pos = strtol(ele, &tail, 0);
    }
    else if (i == 3) {
      info -> seq = (ele[0] >= 'a' && ele[0] <= 'z') ? ele[0] + 32 : ele[0];
    }
    else if (i == 4) {
      char * tail;
      info -> dep[0] = strtol(ele, &tail, 0);
    }
  }
  return;
}

double depratio(Info * start, long n) {
  double * rat = xmalloc(n * sizeof(double));
  long i = 0;
  Info * p = start;
  while (p) {
    if (p -> ratio > 0) {
      rat[i++] = p -> ratio;
    }
    p = p -> next;
  }
  gsl_sort(rat, 1, i);
  double res = (i % 2 == 1) ? rat[(i - 1) / 2] : (rat[(i - 1) / 2] + rat[i / 2]) / 2;
  xfree(rat);
  return log2(res);
}

double gc(Info * start, long n) {
  int ngc = 0;
  int nat = 0;
  Info * p = start;
  while (p) {
    if (p -> ratio > 0) {
      if (p -> seq == 'G' || p -> seq == 'C') {
	ngc++;
      }
      else if (p -> seq == 'A' || p -> seq == 'T') {
	nat++;
      }
    }
    p = p -> next;
  }
  return (ngc || nat) ? (double)ngc / (ngc + nat) : 0;
}

Param * getoption(int argc, char *argv[]){
  Param * param = xmalloc(sizeof(Param));
  param -> w = 1000000;
  param -> s = 1000;
  param -> d = 10;
  param -> n = 100;
  memset(param -> file1, 0, 1000);
  memset(param -> file2, 0, 1000);
  char c;
  extern char *optarg;
  extern int optind, optopt;
  while ((c = getopt(argc, argv, "i:w:s:n:d:f:F:p:c:t:I:h")) != (char)-1) {
    switch(c) {
      case 'w':
        if (sscanf(optarg, "%ld", &param -> w) != 1)
          usage(argv[0]);
        break;
      case 's':
        if (sscanf(optarg, "%ld", &param -> s) != 1)
          usage(argv[0]);
        break;
      case 'd':
        if (sscanf(optarg, "%d", &param -> d) != 1)
          usage(argv[0]);
        break;
      case 'n':
        if (sscanf(optarg, "%d", &param -> n) != 1)
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
           -s : step size [1000]      \n\
           -d : maximum allele frequency [0.95]      \n\
           -n : maximum iteration [10]\n             \n\
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
