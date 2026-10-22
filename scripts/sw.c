/* Smith-Waterman (Gotoh affine gaps, BLOSUM62, gap open 11 / extend 1 = BLAST defaults) : query set vs database.
   Output: query \t target \t raw_score \t qlen \t tlen \t q_end \t t_end
   Usage: sw queries.fa db.fa min_score > out.tsv  */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
static const char *AA="ARNDCQEGHILKMFPSTWYV";
static int B[20][20]={
{4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0},{-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3},
{-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3},{-2,-2,1,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3},
{0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1},{-1,1,0,0,-3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2},
{-1,0,0,2,-4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2},{0,-2,0,-1,-3,-2,-2,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3},
{-2,0,1,-1,-3,0,0,-2,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3},{-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3,-1,3},
{-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1},{-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2},
{-1,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,1},{-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1},
{-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2},{1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2},
{0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0},{-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3},
{-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1},{0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4}};
static int S[21][21]; static int idx[256];
typedef struct { char *name; unsigned char *s; int n; } Seq;
static Seq *readfa(const char *p,int *cnt){
  FILE *f=fopen(p,"r"); if(!f){perror(p);exit(1);} int cap=1024,n=0; Seq *a=malloc(cap*sizeof(Seq));
  char *line=NULL; size_t L=0; ssize_t r; int scap=0; unsigned char *buf=NULL; int bl=0; char *nm=NULL;
  while((r=getline(&line,&L,f))>=0){
    if(line[0]=='>'){
      if(nm){ a[n-1].s=malloc(bl+1); memcpy(a[n-1].s,buf,bl); a[n-1].n=bl; }
      if(n==cap){cap*=2;a=realloc(a,cap*sizeof(Seq));}
      char *e=line+1; while(*e&&!isspace((unsigned char)*e))e++; *e=0; nm=strdup(line+1); a[n].name=nm; n++; bl=0;
    } else { for(char *c=line;*c;c++){ if(isalpha((unsigned char)*c)){ if(bl+1>=scap){scap=scap?scap*2:65536; buf=realloc(buf,scap);} buf[bl++]=(idx[(unsigned char)toupper(*c)]); } } }
  }
  if(nm){ a[n-1].s=malloc(bl+1); memcpy(a[n-1].s,buf,bl); a[n-1].n=bl; }
  fclose(f); *cnt=n; return a; }
int main(int argc,char**argv){
  if(argc<4){fprintf(stderr,"usage\n");return 1;}
  for(int i=0;i<256;i++)idx[i]=20; for(int i=0;i<20;i++)idx[(int)AA[i]]=i;
  for(int i=0;i<21;i++)for(int j=0;j<21;j++)S[i][j]=(i<20&&j<20)?B[i][j]:-1;
  int nq,nd; Seq *Q=readfa(argv[1],&nq),*D=readfa(argv[2],&nd); int minsc=atoi(argv[3]);
  const int GO=11, GE=1; int maxl=0; for(int i=0;i<nd;i++) if(D[i].n>maxl)maxl=D[i].n;
  int *H=malloc((maxl+2)*sizeof(int)), *E=malloc((maxl+2)*sizeof(int));
  for(int q=0;q<nq;q++){
    int m=Q[q].n; int (*prof)[21]=malloc(sizeof(int)*21*(m)); 
    for(int i=0;i<m;i++)for(int a=0;a<21;a++)prof[i][a]=S[Q[q].s[i]][a];
    for(int d=0;d<nd;d++){
      int n=D[d].n; unsigned char *t=D[d].s; int best=0,bi=0,bj=0;
      for(int j=0;j<=n;j++){H[j]=0;E[j]=0;}
      for(int i=1;i<=m;i++){
        int diag=0, F=0, hleft=0; int *pr=prof[i-1];
        for(int j=1;j<=n;j++){
          int e=E[j]-GE; int e2=H[j]-GO-GE; if(e2>e)e=e2; if(e<0)e=0; E[j]=e;
          int f=F-GE; int f2=hleft-GO-GE; if(f2>f)f=f2; if(f<0)f=0; F=f;
          int h=diag+pr[t[j-1]]; if(e>h)h=e; if(f>h)h=f; if(h<0)h=0;
          diag=H[j]; H[j]=h; hleft=h;
          if(h>best){best=h;bi=i;bj=j;}
        }
      }
      if(best>=minsc) printf("%s\t%s\t%d\t%d\t%d\t%d\t%d\n",Q[q].name,D[d].name,best,m,n,bi,bj);
    }
    free(prof);
  }
  return 0; }
