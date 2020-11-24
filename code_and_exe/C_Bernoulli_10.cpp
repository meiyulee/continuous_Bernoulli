#include <math.h>
#include <stdio.h>
#include <stdlib.h> 
#include <malloc.h>
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
 // input the lamda to convert the expected value
double Continuous_Bernoulli_1(double X)
{ 
   double KK[22],tep1,tep2=-1;
   int i;
   if (X>=0.001 && X<=0.999)
   {
   	tep1=0.279390+0.441311*X-0.500045730711711430000000000000;
   	KK[1]=tep1;
   	for(i=2;i<=21;++i) KK[i]=KK[i-1]*tep1;
   	tep2=0.500058872934914690000000000000+
   	     0.773594830650836230000000000000*KK[1]+
   	     -0.015152112930081785000000000000*KK[2]+
   	     -27.279009342193604000000000000000*KK[3]+
   	     10.363707900047302000000000000000*KK[4]+
   	     15822.388427734375000000000000000000*KK[5]+
   	     -2817.424682617187500000000000000000*KK[6]+
   	     -3612752.687500000000000000000000000000*KK[7]+
   	     391281.722656250000000000000000000000*KK[8]+
   	     452401608.000000000000000000000000000000*KK[9]+
   	     -31440996.250000000000000000000000000000*KK[10]+
   	     -33874673664.000000000000000000000000000000*KK[11]+
   	     1540792624.000000000000000000000000000000*KK[12]+
   	     1582581137408.000000000000000000000000000000*KK[13]+
   	     -46642316288.000000000000000000000000000000*KK[14]+
   	     -46495537037312.000000000000000000000000000000*KK[15]+
   	     850124546048.000000000000000000000000000000*KK[16]+
   	     834533872107520.000000000000000000000000000000*KK[17]+
   	     -8542741594112.000000000000000000000000000000*KK[18]+
   	     -8357328558489600.000000000000000000000000000000*KK[19]+
   	     36339642531840.000000000000000000000000000000*KK[20]+
   	     35775834451083264.000000000000000000000000000000*KK[21];
   }
  return(tep2);
 }
 // input the expected value to convert the variance
double Continuous_Bernoulli_2(double X)
{ 
   double KK[22],tep1,tep2=-1;
   int i;
   if (X>=0.143853919 && X<=0.522534496)
   {
   	tep1=0.019991+0.134196*X-0.074925030843033424000000000000;
   	KK[1]=tep1;
   	for(i=2;i<=6;++i) KK[i]=KK[i-1]*tep1;
   	tep2= 0.078453112965998995000000000000+
          0.798566816470611230000000000000*KK[1]+
          -31.711886927382693000000000000000*KK[2]+
          -96.513433280750178000000000000000*KK[3]+
          2025.523392695933600000000000000000*KK[4]+
          -5401.664180755615200000000000000000*KK[5]+
          1123298.764297485400000000000000000000*KK[6];
   }
   if (X>=0.522534496 && X<=0.856221427)
   {
   	tep1=0.163707+-0.148751*X-0.073800792252000508000000000000;
   	KK[1]=tep1;
   	for(i=2;i<=14;++i) KK[i]=KK[i-1]*tep1;
   	tep2= 0.076871792010223317000000000000+
          0.826624807828920890000000000000*KK[1]+
          -25.39717203173351000000000000000*KK[2]+
          -157.130625724792480000000000000000*KK[3]+
          2458.807250976562500000000000000000*KK[4]+
          1185575.781250000000000000000000000000*KK[5]+
          4707383.000000000000000000000000000000*KK[6]+
          -8267723008.000000000000000000000000000000*KK[7]+
          -185558838784.000000000000000000000000000000*KK[8]+
          25159116259328.000000000000000000000000000000*KK[9]+
          1107202169372672.000000000000000000000000000000*KK[10]+
          -16725473455243264.000000000000000000000000000000*KK[11]+
          -1902815785101819900.000000000000000000000000000000*KK[12]+
          -41867494033524261000.000000000000000000000000000000*KK[13]+
          -307130787519134700000.000000000000000000000000000000*KK[14];
   }
  return(tep2);
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

// Z distribution database and the size is n
void Zdatabase(int n,double *datastore)
{
	int i,j;
	double tep;
	for(i=0;i<n;++i)
	{
	   tep=simulated_normal(0,1);
	   datastore[i]=tep;
	}
	divide_Conquer_X_step(datastore,n);
}
// the test statistic sampling distribution when the sample sizes are samplesize1 and samplesize2
double *test_sampling_2_t_spool(int samplesize1,double *sampledata1,int samplesize2,double *sampledata2,double lamda)
{
	double *datastore,samplemean,samplemean1,samplemean2,VarX,sigma;
	int i,j,n=1000000,samplesize;	
	
	datastore=(double *)malloc(n*sizeof(double));
	samplesize=samplesize1+samplesize2;
    for(i=0;i<n;++i)
    { 
    	samplemean1=0;
    	for(j=0;j<samplesize1;++j)
    	{
    		sampledata1[j]=simulated_continue_bernoulli(lamda);
    		samplemean1=samplemean1+sampledata1[j];
		}
	   samplemean=samplemean1;
       samplemean1=samplemean1/(double)samplesize1;
       samplemean2=0;
    	for(j=0;j<samplesize2;++j)
    	{
    		sampledata2[j]=simulated_continue_bernoulli(lamda);
    		samplemean2=samplemean2+sampledata2[j];
		}
	   samplemean=samplemean+samplemean2;
       samplemean2=samplemean2/(double)samplesize2;
       samplemean=samplemean/(double)samplesize;
       VarX=0;
       for(j=0;j<samplesize1;++j)
       VarX=VarX+(sampledata1[j]-samplemean)*(sampledata1[j]-samplemean);
       for(j=0;j<samplesize2;++j)
       VarX=VarX+(sampledata2[j]-samplemean)*(sampledata2[j]-samplemean);
       VarX=VarX/(double)(samplesize-1); 
       sigma=exp(0.5*log(VarX/(double)samplesize1+VarX/(double)samplesize2));
       datastore[i]=(samplemean1-samplemean2)/sigma;
	}
	divide_Conquer_X_step(datastore,n);
	return(datastore);
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
     double *datastoreX,*sampledata1,*sampledata2,samplemean1,samplemean2,samplemean;
     double pvalue,EXX,VarXX,VarX,sigma,c1,c2,lamda,lamda1,lamda2,VarX1,VarX2,Ztest;
	 int i,samplesize=0,samplesize1=1,samplesize2=1,j=0,Zflag=0;
	 char *ss,*filename;
	 FILE *f;
	
      ss=(char *)malloc(1024*sizeof(char));
      filename=(char *)malloc(1024*sizeof(char));
      loop1:
      printf(" Please input the 1st data filename>");
      scanf("%s",filename);
      printf("%s\n",filename);
	  f=fopen(filename,"r");
	  if (f==NULL) 
	  { 
	   printf(" %s cannot be found, please input filename again.\n",filename);
	   goto loop1;
	  }
	  
	  loop2:
	  j=get_string(f,ss);
	  if (j>=0) { samplesize1++; goto loop2;}	
	  if (j==-1) samplesize1--;
	  fclose(f);
	  
	  printf(" 1st sample data size=%d\n",samplesize1);
	  f=fopen(filename,"r");
      sampledata1=(double *)malloc(samplesize1*sizeof(double));
      for(j=0;j<samplesize1;++j) 
       {  get_string(f,ss);
          sampledata1[j]=atof(ss);
          if (sampledata1[j]<0 || sampledata1[j]>1) 
          {
          	printf(" The %d-th sample data of 1st file is error,\n which value=%f is not in [0,1],\n this program will be terminated.",j+1,sampledata1[j]);
		    goto end;
		  }
	   }
	  fclose(f);
      
      loop3:
      printf(" Please input the 2nd data filename>");
      scanf("%s",filename);
      printf("%s\n",filename);
	  f=fopen(filename,"r");
	  if (f==NULL) 
	  { 
	   printf(" %s cannot be found, please input filename again.\n",filename);
	   goto loop3;
	  }
	  
	  loop4:
	  j=get_string(f,ss);
	  if (j>=0) { samplesize2++; goto loop4;}	
	  if (j==-1) samplesize2--;
	  fclose(f);
	  
	  printf(" 2nd sample data size=%d\n",samplesize2);
	  f=fopen(filename,"r");
      sampledata2=(double *)malloc(samplesize2*sizeof(double));
      for(j=0;j<samplesize2;++j) 
       {  get_string(f,ss);
          sampledata2[j]=atof(ss);
          if (sampledata2[j]<0 || sampledata2[j]>1) 
          {
          	printf(" The %d-th sample data of 2nd file is error,\n which value=%f is not in [0,1],\n this program will be terminated.",j+1,sampledata2[j]);
		    goto end;
		  }
	   }
	  fclose(f);
	  samplemean=0; samplesize=samplesize1+samplesize2;
      for(i=0;i<samplesize1;++i)
        samplemean=samplemean+sampledata1[i];
      samplemean1=samplemean/(double)(samplesize1);
      samplemean2=0;
      for(i=0;i<samplesize2;++i)
      {
      	samplemean=samplemean+sampledata2[i];
      	samplemean2=samplemean2+sampledata2[i];
	  }
	  samplemean2=samplemean2/(double)(samplesize2);
	  samplemean=samplemean/(double)(samplesize);
	  lamda=Continuous_Bernoulli_pointer(samplemean);
	  lamda1=Continuous_Bernoulli_pointer(samplemean1);
	  lamda2=Continuous_Bernoulli_pointer(samplemean2);
	  if (lamda1<0 || lamda2<0 || lamda<0) { printf(" H:lamda1=lamda2 cannot be testing\n"); goto end;}
	  samplemean2=0;
      for(i=0;i<samplesize2;++i)
        samplemean2=samplemean2+sampledata2[i];
	  samplemean2=samplemean2/(double)(samplesize2);
	  VarX1=0;
	  for(i=0;i<samplesize1;++i)
       VarX1=VarX1+(sampledata1[i]-samplemean1)*(sampledata1[i]-samplemean1);
      VarX1=VarX1/(double)(samplesize1-1);
      VarX2=0;
	  for(i=0;i<samplesize2;++i)
       VarX2=VarX2+(sampledata2[i]-samplemean2)*(sampledata2[i]-samplemean2);
      VarX2=VarX2/(double)(samplesize2-1);
      
	  VarX=0;
      for(i=0;i<samplesize1;++i)
       VarX=VarX+(sampledata1[i]-samplemean)*(sampledata1[i]-samplemean);
      for(i=0;i<samplesize2;++i)
       VarX=VarX+(sampledata2[i]-samplemean)*(sampledata2[i]-samplemean);
      VarX=VarX/(double)(samplesize-1); 
      sigma=exp(0.5*log(VarX/(double)samplesize+VarX/(double)samplesize));
      printf(" The 1st data,......\n");
	  printf(" sample size=%d, sample mean=%20.10f\n",samplesize1,samplemean1);
      printf(" sample variance=%20.10f,\n lamda1 estimated value=%20.10f\n",VarX1,lamda1);
      printf(" The 2nd data,......\n");
	  printf(" sample size=%d, sample mean=%20.10f\n",samplesize2,samplemean2);
      printf(" sample variance=%20.10f,\n lamda2 estimated value=%20.10f\n",VarX2,lamda2);
      printf(" Under H0:lamda1=lamda2,.....\n");
      printf(" sample size=%d, sample mean=%20.10f\n",samplesize,samplemean);
      printf(" sample variance=%20.10f,\n lamda1 estimated value=%20.10f\n",VarX,lamda);
	  if (0.1<=lamda &&  lamda<=0.9)
      {
      	if (lamda<0.5)
      	i=33+(int)((double)350*(0.5-lamda));
      	else 
      	i=33+(int)((double)350*(lamda-0.5));
	  }
	  if (0.1>lamda) i=500+(int)((double)15000*(0.1-lamda));
	  if (lamda>0.9) i=500+(int)((double)15000*(lamda-0.9));  
	  if (samplesize>=i) Zflag=1;
	  else Zflag=0;
	  printf(" sample size=%d, the requirement size=%d for Z distribution,\n",samplesize,i);
      
      switch (Zflag)
	  {
	  case 0:
	  printf(" please wait a moment to collect the sampling distribution\n");
	  datastoreX=test_sampling_2_t_spool(samplesize1,sampledata1,samplesize2,sampledata2,lamda);
	  break;
      case 1:
      datastoreX=(double *)malloc(1000000*sizeof(double));
      Zdatabase(1000000,datastoreX);
      break;
      }
      printf(" H0:lamda1=lamda2,\n");
      EXX=Continuous_Bernoulli_1(lamda);
      VarXX=Continuous_Bernoulli_2(EXX);
      printf(" The population estimated mean=%20.10f,\n population estimated variance=%20.10f under H0\n",EXX,VarXX);
	  if (lamda>=0)
	  {
      Ztest=(samplemean1-samplemean2)/sigma;
      for(i=0;i<1000000;++i) if (Ztest<datastoreX[i]) break;
      pvalue=(double)(i+1)/(double)1000000;
      if (pvalue>=1) pvalue=1;
      if (pvalue<0.5) pvalue=2*pvalue; else pvalue=2*(1-pvalue);
      printf(" test statistic=%20.10f, p value=%f,\n",Ztest,pvalue);
      if (pvalue>=0.05) // the H0:mu1=mu2c is not rejected, estimated lamda1 and lamda2
        printf(" lamda1 hat=lamda2 hat=lamda hat=%f hen H0:lamda1=lamda2 not rejected.\n",lamda);
      else
        printf(" lamda1 hat=%f, lamda2 hat=%f hen H0:lamda1=lamda2 rejected.\n",lamda1,lamda2);
    	  }
      else printf(" The test cannot be done when lamda hat cannot be found.\n");
      end:
	  system("pause");
	  free(filename); free(ss); free(datastoreX);
	  return 0;
}
