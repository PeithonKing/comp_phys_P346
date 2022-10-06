// Cholesky decomposition

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int is_symmetric(float **a,int n);
void cholesky(float **a,int n);
void forward(float **a,float *b,float *y,int n);
void backward(float **a,float *y,float *x,int n);

int main(int argc,char *argv[])
{

  if(argc!=4){
    printf("Usage : ./a.out <matrix-dim> <in-matrix-file> <out-matrix-file>\n"); exit(1);
  }

  // open files
  FILE *fin,*fout;
  fin=fopen(argv[2],"r");
  fout=fopen(argv[3],"w");

  // allocate space
  int i,j,N;
  float **A,*B;
  float *Y,*X;
  N=atoi(argv[1]);
  A=(float **)malloc(N*sizeof(float *));
  for(i=0;i<N;i++) A[i]=(float *)malloc(N*sizeof(float));
  B=(float *)malloc(N*sizeof(float));
  Y=(float *)malloc(N*sizeof(float));
  X=(float *)malloc(N*sizeof(float));

  // read data
  int ctr=0;
  float temp; int offset;
  char strng[512];
  while(fgets(strng,sizeof(strng),fin)!=NULL){
     if(strng[0]=='#'){ i=j=0; continue; }
     char *data=strng;
     if(ctr<N){while(sscanf(data,"%f%n",&temp,&offset)==1){
	A[i][j]=temp; j++; data+=offset;
     }}
     else{while(sscanf(data,"%f%n",&temp,&offset)==1){
        B[i]=temp; data+=offset;
     }}
     i++;j=0;ctr++;
  } // end of while

  // print read-in matrix
  fprintf(fout,"# Original matrix\n");
  for(i=0;i<N;i++){
     for(j=0;j<N;j++) fprintf(fout,"%.2f  ",A[i][j]);
     fprintf(fout,"\n");
  }
  fprintf(fout,"\n# RHS\n");
  for(i=0;i<N;i++) fprintf(fout,"%.2f\n",B[i]);

  // check data and cholesky decompose
  int check;
  check=is_symmetric(A,N);
  if(check==1){ fprintf(fout,"Matrix not symmetric! Exiting!!\n"); exit(1); }
  else{ cholesky(A,N); }

  // print cholesky matrix
  fprintf(fout,"\nCholesky decomposed matrix!\n");
  for(i=0;i<N;i++){
     for(j=0;j<N;j++){
        fprintf(fout,"%.2f  ",A[i][j]);
     }
     fprintf(fout,"\n");
  }

  // forward - backward substitution
  for(i=0;i<N;i++) Y[i]=X[i]=0.0;
  forward(A,B,Y,N);
  fprintf(fout,"\nForward solution Y:\n");
  for(i=0;i<N;i++) fprintf(fout,"%.f\n",Y[i]);
  backward(A,Y,X,N);

  fprintf(fout,"\nSolution is :\n");
  for(i=0;i<N;i++) fprintf(fout,"%.f\n",X[i]);


} // end of main

int is_symmetric(float **a,int n)
{

  int l,m;

  for(l=0;l<n;l++)for(m=l+1;m<n;m++){
     if(a[l][m]!=a[m][l]) return(1);
  }

  return(2);

} // end of is_symmetric()

void cholesky(float **a,int n)
{

  int i,j,k;
  float temp;

  for(i=0;i<n;i++){

     temp=0;
     for(j=0;j<i;j++) temp+=a[j][i]*a[j][i];
     a[i][i]=(float)sqrt(a[i][i]-temp);

     for(j=i+1;j<n;j++){
	temp=0;
	for(k=0;k<i;k++) temp+=a[i][k]*a[k][j];
	a[j][i]=(a[j][i]-temp)/a[i][i];
	a[i][j]=a[j][i];
	}

     } // end of i-loop

} // end of cholesky()

void forward(float **a,float *b,float *y,int n)
{

  int i,j;
  float sum;

  for(i=0;i<n;i++){
     sum=0;
     for(j=0;j<i;j++) sum+=a[i][j]*y[j];
     y[i]=(b[i]-sum)/a[i][i];
     }

} // end of forward()

void backward(float **a,float *y,float *x,int n)
{

  int i,j;
  float sum;

  for(i=n-1;i>=0;i--){
     sum=0;
     for(j=i+1;j<n;j++) sum+=a[i][j]*x[j];
     x[i]=(y[i]-sum)/a[i][i];
     }

} // end of forward()

