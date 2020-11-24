#include <math.h>
#include <stdio.h>
#include <stdlib.h> 
#include <malloc.h>
// input the smaple mean to convert the lamda point estimated value
double Continuous_Bernoulli_pointer(double X)
{ 
   double KK[20],tep1,tep2=-1; 
// if the sample mean is not in [0.143853919,0.856221427] return -1, it is error
   int i;
   if (X>=0.143853919 && X<=0.856221427)
   {
   	tep1=-0.596698+2.193196*X-0.500000000000000560000000000000;
   	KK[1]=tep1;
   	for(i=2;i<=18;++i) KK[i]=KK[i-1]*tep1;
   	tep2=0.499973865804236080000000000000+
         1.368024096854642700000000000000*KK[1]+
         -0.000924747670069336890000000000*KK[2]+
         -2.736078237077606400000000000000*KK[3]+
         0.095109043642878532000000000000*KK[4]+
         5.748377367592183900000000000000*KK[5]+
         -1.841998845338821400000000000000*KK[6]+
         -12.357242575206328000000000000000*KK[7]+
         16.361405849456787000000000000000*KK[8]+
         26.417928500100970000000000000000*KK[9]+
         -80.021261215209961000000000000000*KK[10]+
         -48.621550429612398000000000000000*KK[11]+
         228.768722534179690000000000000000*KK[12]+
         64.702439151704311000000000000000*KK[13]+
         -380.758743286132810000000000000000*KK[14]+
         -51.895506033673882000000000000000*KK[15]+
         341.663604736328120000000000000000*KK[16]+
         18.360968290828168000000000000000*KK[17]+
         -127.708103179931640000000000000000*KK[18];
   }
  return(tep2);
 }
 // input the smaple data and computing the sample mean
 // the sample mean to convert the lamda point estimated value
double Continuous_Bernoulli_pointer(double *sample,int samplesize,double *lamda)
{
	double sumX;
	int i,j=0;
	
	sumX=0; 
    for(i=0;i<samplesize;++i)
       	sumX=sumX+sample[i];
	sumX=sumX/(double)samplesize;
	again:
	*lamda=Continuous_Bernoulli_pointer(sumX);
	if (*lamda<0.001) *lamda=0;
	if (*lamda>=0.999) *lamda=1;
	return(sumX);
}
// sorting data
int divide_Conquer_X_step(double *doubledata, int datasize)
{
  int i,count1,count2;
  double MAX,MIN,compare_point;
  double *tepdata1,*tepdata2;
  
  if (datasize<=1) return(0); 
  MAX=doubledata[0];
  MIN=doubledata[0];
  
  for(i=1;i<datasize;++i) 
  {
  	if (doubledata[i]>MAX) MAX=doubledata[i]; 
	if (doubledata[i]<MIN) MIN=doubledata[i]; 
  }
  compare_point=(MIN+MAX)/(double)2;
  tepdata1=(double *)malloc(datasize*sizeof(double));
  tepdata2=(double *)malloc(datasize*sizeof(double));
  count1=0;  count2=0;    
  for(i=0;i<datasize;++i) 
  {
  	if (doubledata[i]<compare_point) 
  	{ 
  	  tepdata1[count1]=doubledata[i];
  	  ++count1;
	  }
	else  
	{ 
  	  tepdata2[count2]=doubledata[i];
  	  ++count2;
	  }
  }
  if (count1>=1) divide_Conquer_X_step(tepdata1,count1);
  if (count2>=1) divide_Conquer_X_step(tepdata2,count2);
 
  for(i=0;i<count1;++i)
    doubledata[i]=tepdata1[i];
  for(i=0;i<count2;++i)
    doubledata[i+count1]=tepdata2[i];
  free(tepdata1);  
  free(tepdata2);   
}
// the random number 
double simulated_uniform(double alpha,double beta)
{
   double tep1,tep=0,dividewangx;
 
   dividewangx=(double)32767*(double)32767*(double)32767;
   again:
   tep1=((double)rand()*(double)32767*(double)32767+(double)rand()*(double)32767+(double)rand())/dividewangx;
   if (tep1<0 || tep1>1) goto again;
   tep=alpha+tep1*(beta-alpha);   
   return(tep);
}
// the normal distribution simulator
double simulated_normal(double mu,double sigma)
{
   double tep1,tep2,tep=0,tepx,pii=M_PI;

   again02:
   tep1=simulated_uniform(0,1); 
   tep2=simulated_uniform(0,1);
   tepx=-2*log(tep1);
   if (tepx>0) tepx=exp(0.5*log(tepx)); else tepx=0;
   if (rand()<16384) tep=tepx*sin(2*tep2*pii); 
         else tep=tepx*cos(2*tep2*pii); 
   tep=mu+tep*sigma;
   return(tep);
}
// the chi square distribution simulator
void chi_square_simulator(int df,int n,double *datastore)
{
	int i,j;
	double tep,sum;
	for(i=0;i<n;++i)
	{
	   sum=0;
	   for(j=0;j<df;++j)
	   {
	   	tep=simulated_normal(0,1);
	   	sum=sum+tep*tep;
	   }
	   datastore[i]=sum;
	}
	divide_Conquer_X_step(datastore,n);
}
// The comparison of two strings and the index upper limit is len
int strcompare_spc(char *s,char *ss,int len)
{
  int i=0;
  for(;s[i]!=0 && ss[i]!=0 && i<len;++i)
   if (s[i]!=ss[i]) break;
  return(s[i]-ss[i]);
// return 0 when two strings are eqully
}
// Continuous Bernoulli distribution simulator
double simulated_continue_bernoulli(double lamda)
{
   double tep=0,tep1;
   char ss[1024];
   int kk=1;
    
    if (lamda>=0.499 && lamda<=0.501)
    {
     sprintf(ss,"%12.10f",lamda);
     kk=strcompare_spc(ss,(char *)"0.5000000000",12);
	}
	if (kk!=0) 
    {
         tep1=simulated_uniform(0,1);
         tep1=(2*lamda-1)*tep1-lamda+1;
         tep1=tep1/(1-lamda);
         tep=log(tep1)/log(lamda/(1-lamda));
		}
		else
		  tep=simulated_uniform(0,1);
   return(tep);
}

// Continuous Bernoulli distribution simulator giving the cumulative probability
double ret_continuous_bernoulli_x(double lamda,double pr)
{
	double tep=0,tep1;
    char ss[1024];
    int kk=0;
    
    sprintf(ss,"%12.10f",lamda);
    kk=strcompare_spc(ss,(char *)"0.5000000000",12);    
    // kk is 0 when lamda=0.5, X~Uniform(0,1) if lamda=0.5
    if (kk!=0) 
    {
         tep1=pr;
         tep1=(2*lamda-1)*tep1-lamda+1;
         tep1=tep1/(1-lamda);
         tep=log(tep1)/log(lamda/(1-lamda));
		}
		else
		  tep=pr;
   return(tep);
}
// degree of freedom
int ret_df(int n)
{
  double tepw1;
  double pvalue;
  int classx;
 
  tepw1=1;
  for (classx=0;;++classx) 
  {
  	tepw1=tepw1*2;
  	if (tepw1>(double)n) break;
  }
  classx++;  
  return(classx-2);
}
// goodness of fit
double goodness_of_fitX(int n,double *data,int *dfret,double lamda,double *ret_pvalue,double *datastore)
{
  double max,min,*upper,*lower,*expected,*actual,tepw,tepw1,pr,prcum=0,*getvalue1,*teptep;
  double pvalue;
  int classx,i,j,flagx=0,flagx1=0,nxn=10000000;
  char *s;
  FILE *f;
  
  tepw1=1;
  for (classx=0;;++classx) 
  {
  	tepw1=tepw1*2;
  	if (tepw1>(double)n) break;
  }
  classx++;  
  upper=(double *)malloc(2000*sizeof(double));
  lower=(double *)malloc(2000*sizeof(double));
  expected=(double *)malloc(2000*sizeof(double));
  actual=(double *)malloc(2000*sizeof(double));
  pr=1/(double)classx;
  if (pr<0.05 && pr*(double)n<5) {pr=0.05; classx=20;  } 
  *dfret=classx-1; 
  for(i=0;i<classx;++i) {expected[i]=(double)n/(double)classx;  actual[i]=0;}
  lower[0]=0;
  for(prcum=pr,i=1;i<=classx;prcum=prcum+pr,++i)
  {
   upper[i-1]=ret_continuous_bernoulli_x(lamda,prcum); lower[i]=upper[i-1];}
   upper[classx-1]=1;
   for(j=0;j<n;++j) 
  { for(i=0;i<classx;++i)
    { if (i==0 && data[j]<=upper[i]) actual[i]=actual[i]+1;
      if (i>0 && data[j]<=upper[i]&& data[j]>lower[i]) actual[i]=actual[i]+1;
    }
  }
  actual[classx-1]=n;
  for(i=0;i<classx-1;++i) actual[classx-1]=actual[classx-1]-actual[i];
  tepw=0;
  for(i=0;i<classx;++i) tepw=tepw+(expected[i]-actual[i])*(expected[i]-actual[i])/expected[i]; 
  for(i=0;i<1000000;++i) if (datastore[i]>=tepw) break;
  pvalue=(double)(i+1)/(double)1000000;
  if (pvalue<0) pvalue=0;
  if (pvalue>1) pvalue=1;
  printf("\n  pearson goodness of fit\n");
  printf(" H0:X~Continuous Bernoulli(lamda),\n");
  printf(" H1:against H0,\n");
  printf(" lamda under H0=%20.10f\n",lamda);
  printf("  pearson goodness of fit\n");
  printf("\n degree of freedom=%d\n",*dfret);
  printf(" chi square test=%f\n",tepw);
  printf(" p value=%f\n",1-pvalue);
  *ret_pvalue=1-pvalue;
  end:
  free(upper); free(lower); free(expected); free(actual);
  return(tepw);
}
// degree of freedom
int ret_dfx(int n)
{
  double tepw1;
  double pvalue;
  int classx;
 
  tepw1=1;
  for (classx=0;;++classx) 
  {
  	tepw1=tepw1*2;
  	if (tepw1>(double)n) break;
  }
  classx++;  
  return(classx-1);
}
// reading the string from file
int get_string(FILE *f,char *s)
{
  int i,ret=0;
  char c;
  for(i=0;;++i)
  {
    c=fgetc(f);
    if (c!='\n' && c!=EOF )  s[i]=c;
    else break;
   }
   s[i]=0;
   if (c==EOF) ret=-1+-1*i;
   // -1 is EOF and no data
   // less than is EOF and has data
   return(ret);
}
int main() {  

     double *sampledata,tep,lamdahat,samplemean,pvalue,*datastoreX,p1;
	 int i,j,samplesize=1,df;
	 char *filename,*ss;
	 FILE *f;

      datastoreX=(double *)malloc(1000000*sizeof(double));
      ss=(char *)malloc(1024*sizeof(char));
      filename=(char *)malloc(1024*sizeof(char));
      // input the filename which has the sample data
      loop1:
      printf(" Please input the data filename:");
      scanf("%s",filename);
      printf("%s\n",filename);
	  f=fopen(filename,"r");
	  if (f==NULL) 
	  { 
	   printf(" %s cannot be found, please input filename again.\n",filename);
	   goto loop1;
	  }
	  // computing sample data size
	  loop2:
	  j=get_string(f,ss);
	  if (j>=0) { samplesize++; goto loop2;}	
	  if (j==-1) samplesize--;
	  fclose(f);
	   // reading sample data
	  printf(" sample data size=%d\n",samplesize);
	  f=fopen(filename,"r");
      sampledata=(double *)malloc(samplesize*sizeof(double));
      for(j=0;j<samplesize;++j) 
       {  get_string(f,ss); 
          sampledata[j]=atof(ss);
          if (sampledata[j]<0 || sampledata[j]>1) 
          {
          	printf(" The %d-th sample data is error,\n which value=%f is not in [0,1],\n this program will be terminated.",j+1,sampledata[j]);
		    goto end;
		  }
	   }
	  fclose(f);
	  
	  // input lamda value
	  loop3:
	  printf(" H0:Continuous Bernoulli(lamda),lamda=");
	  scanf("%s",ss);
	  p1=atof(ss);
	  if (p1<=0.001 || p1>=0.999)
	  {
	  	printf(" The lamda=%f has a problem, which range is [0.001,0.999],\n please input lamda again.",p1);
	  	goto loop3;
	  }
      df=ret_dfx(samplesize); // degree of freedom when doing goodness of fit
      printf(" Collecting chi square(df=%d) database,...\n",df);
	  chi_square_simulator(df,1000000,datastoreX);
      // computing sample mean and estimated lamda 
	  samplemean=Continuous_Bernoulli_pointer(sampledata,samplesize,&lamdahat);
	  printf(" the sample mean=%20.10f\n the lamda estimated=%f.\n",samplemean,lamdahat);
	  if (lamdahat<0) printf(" the sample mean=%20.10f\n the lamda hat is error.\n",samplemean);
	  // the goodness of fit
      goodness_of_fitX(samplesize,sampledata,&df,p1,&pvalue,datastoreX);
      end:
	  system("pause");
      free(sampledata); free(filename);
	  free(ss);  free(datastoreX);
	  return 0;
}
      
      
