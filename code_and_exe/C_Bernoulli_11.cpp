#include <math.h>
#include <stdio.h>
#include <stdlib.h> 
#include <malloc.h>
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
    
    if (lamda>=0.49 && lamda<=0.51)
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


double Continuous_Bernoulli_pointer(double X)
{ 
   double KK[20],tep1,tep2=-1;
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

double *test_sampling_2_t(int samplesize,double *sampledata,int samplesize2,double *sampledata2,double lamda,double lamdaY)
{
	double *datastore,samplemean,EX,VarX,samplemean2,EX2,VarX2,sigmaXbar;
	int i,j,n=1000000;	
	
	EX=Continuous_Bernoulli_1(lamda);
	EX2=Continuous_Bernoulli_1(lamdaY);
	datastore=(double *)malloc(n*sizeof(double));
    for(i=0;i<n;++i)
    {
    	samplemean=0;VarX=0;
    	for(j=0;j<samplesize;++j)
    	{
    		sampledata[j]=simulated_continue_bernoulli(lamda);
    		samplemean=samplemean+sampledata[j];
		}
       samplemean=samplemean/(double)samplesize;
       VarX=0;
       for(j=0;j<samplesize;++j)
       VarX=VarX+(sampledata[j]-samplemean)*(sampledata[j]-samplemean);
       VarX=VarX/(double)(samplesize-1); 
       samplemean2=0;VarX2=0;
    	for(j=0;j<samplesize2;++j)
    	{
    		sampledata2[j]=simulated_continue_bernoulli(lamdaY);
    		samplemean2=samplemean2+sampledata2[j];
		}
       samplemean2=samplemean2/(double)samplesize2;
       VarX2=0;
       for(j=0;j<samplesize2;++j)
       VarX2=VarX2+(sampledata2[j]-samplemean2)*(sampledata2[j]-samplemean2);
       VarX2=VarX2/(double)(samplesize2-1); 
   	   sigmaXbar=exp(0.5*log(VarX/(double)samplesize+VarX2/(double)samplesize2));
       datastore[i]=(samplemean-samplemean2-(EX-EX2))/sigmaXbar;
	}
	divide_Conquer_X_step(datastore,n);
	return(datastore);
}
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
   return(ret);
}
int main() {  
     double *datastoreX,*sampledata,samplemean,*sampledata2,samplemean2;
     double pvalue,EX,VarX,sigma,EX2,VarX2,c1,c2,lamdaL1,lamdaU1,lamda2,lamda,lamdaY;
	 int i,samplesize=1,samplesize2=1,j=0,Zflag=0,Zflag2=0;
	 char *ss,*filename;
	 FILE *f;
	 double *Xc1,*Xc2;

      Xc1=(double *)malloc(3*sizeof(double));
      Xc2=(double *)malloc(3*sizeof(double));
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
	  if (j>=0) { samplesize++; goto loop2;}	
	  if (j==-1) samplesize--;
	  fclose(f);
	  
	  printf(" 1st sample data size=%d\n",samplesize);
	  f=fopen(filename,"r");
      sampledata=(double *)malloc(samplesize*sizeof(double));
      for(j=0;j<samplesize;++j) 
       {  get_string(f,ss);
          sampledata[j]=atof(ss);
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
	   }
	  fclose(f);
	  samplemean=0;
      for(i=0;i<samplesize;++i)
        samplemean=samplemean+sampledata[i];
	  samplemean=samplemean/(double)(samplesize);
	  lamda=Continuous_Bernoulli_pointer(samplemean);
	  samplemean2=0;
      for(i=0;i<samplesize2;++i)
        samplemean2=samplemean2+sampledata2[i];
	  samplemean2=samplemean2/(double)(samplesize2);
	  lamdaY=Continuous_Bernoulli_pointer(samplemean2);
	  VarX=0;
      for(i=0;i<samplesize;++i)
       VarX=VarX+(sampledata[i]-samplemean)*(sampledata[i]-samplemean);
      VarX=VarX/(double)(samplesize-1); 
      VarX2=0;
      for(i=0;i<samplesize2;++i)
       VarX2=VarX2+(sampledata2[i]-samplemean2)*(sampledata2[i]-samplemean2);
      VarX2=VarX2/(double)(samplesize2-1); 
      sigma=exp(0.5*log(VarX/(double)samplesize+VarX2/(double)samplesize2));
      printf(" The 1st data,......\n");
	  printf(" sample size=%d, sample mean=%20.10f\n",samplesize,samplemean);
	  if (lamda<0) printf(" sample variance=%20.10f,\n lamda estimated value=------\n",VarX);
	  else
      printf(" sample variance=%20.10f,\n lamda estimated value=%20.10f\n",VarX,lamda);
      printf(" The 2nd data,......\n");
	  printf(" sample size=%d, sample mean=%20.10f\n",samplesize2,samplemean2);
	  if (lamdaY<0) printf(" sample variance=%20.10f,\n lamda estimated value=-----\n",VarX2);
	  else
      printf(" sample variance=%20.10f,\n lamda estimated value=%20.10f\n",VarX2,lamdaY);
	  if (lamda<-0.5 || lamdaY<-0.5) { printf(" The C.I. cannot be found\n"); goto end;}
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
	  printf(" 1st sample size=%d, the requirement size=%d for Z distribution,\n",samplesize,i);
	  if (0.1<=lamdaY &&  lamdaY<=0.9)
      {
      	if (lamdaY<0.5)
      	i=33+(int)((double)350*(0.5-lamdaY));
      	else 
      	i=33+(int)((double)350*(lamdaY-0.5));
	  }
	  if (0.1>lamdaY) i=500+(int)((double)15000*(0.1-lamdaY));
	  if (lamdaY>0.9) i=500+(int)((double)15000*(lamdaY-0.9));  
	  if (samplesize2>=i) Zflag2=1;
	  else Zflag2=0;
	  printf(" 2nd sample size=%d, the requirement size=%d for Z distribution,\n",samplesize2,i);
	  if (Zflag==0 || Zflag2==0) 
	  {
	  	Zflag=0;
	  	printf(" the small sample, the testing sampling distribution,\n");
	  }
      else { Zflag=1; 
             printf(" the big sample, the Z distribution,\n");  }
      
      switch (Zflag)
	  {
	  case 0:
	  printf(" please wait a moment to collect the sampling distribution\n");
	  datastoreX=test_sampling_2_t(samplesize,sampledata,samplesize2,sampledata2,lamda,lamdaY);
	  Xc1[0]=datastoreX[(int)((double)1000000*0.05)];
	  Xc2[0]=datastoreX[(int)((double)1000000*0.95)];	
	  Xc1[1]=datastoreX[(int)((double)1000000*0.025)];
	  Xc2[1]=datastoreX[(int)((double)1000000*0.975)];	
	  Xc1[2]=datastoreX[(int)((double)1000000*0.005)];
	  Xc2[2]=datastoreX[(int)((double)1000000*0.995)];	
	  free(datastoreX);
	  break;
      case 1:
      Xc1[0]=-1.645;
      Xc2[0]=1.645;
      Xc1[1]=-1.96;
      Xc2[1]=1.96;
      Xc1[2]=-2.576;
      Xc2[2]=2.576;
      break;
      }
      c1=samplemean-samplemean2+Xc1[0]*sigma;
      c2=samplemean-samplemean+Xc2[0]*sigma;
      printf(" 90%% C.I. for mu1-mu2\n");
      printf(" %20.10f<= mu1- mu2 <= %20.10f\n",c1,c2);
      c1=samplemean-samplemean2+Xc1[1]*sigma;
      c2=samplemean-samplemean2+Xc2[1]*sigma;
      printf(" 95%% C.I. for mu1-mu2\n");
      printf(" %20.10f<= mu1- mu2 <= %20.10f\n",c1,c2);
      c1=samplemean-samplemean2+Xc1[2]*sigma;
      c2=samplemean-samplemean2+Xc2[2]*sigma;
      printf(" 99%% C.I. for mu1-mu2\n");
      printf(" %20.10f<= mu1- mu2 <= %20.10f\n",c1,c2);
      
      c1=samplemean-samplemean2+Xc1[0]*sigma;
      c2=samplemean-samplemean+Xc2[0]*sigma;
      printf(" 90%% C.I. for lamda1-lamda2\n");
      lamda2=Continuous_Bernoulli_pointer(samplemean2);
      lamdaL1=Continuous_Bernoulli_pointer(samplemean+Xc1[0]*sigma);
      lamdaU1=Continuous_Bernoulli_pointer(samplemean+Xc2[0]*sigma);
      if (lamda2<0 || lamdaL1<0 || lamdaU1<0)
	   { printf(" The  90%% C.I. for lamda1-lamda2 cannot be done\n, this program wiil be terminated.");
	     goto end;}
	  printf(" %20.10f<=  lamda1-lamda2  <= %20.10f\n",lamdaL1-lamda2,lamdaU1-lamda2);
	  
	  c1=samplemean-samplemean2+Xc1[1]*sigma;
      c2=samplemean-samplemean+Xc2[1]*sigma;
      printf(" 95%% C.I. for lamda1-lamda2\n");
      lamda2=Continuous_Bernoulli_pointer(samplemean2);
      lamdaL1=Continuous_Bernoulli_pointer(samplemean+Xc1[1]*sigma);
      lamdaU1=Continuous_Bernoulli_pointer(samplemean+Xc2[1]*sigma);
      if (lamda2<0 || lamdaL1<0 || lamdaU1<0)
	   { printf(" The  95%% C.I. for lamda1-lamda2 cannot be done\n, this program wiil be terminated.");
	     goto end;}
	  printf(" %20.10f<=  lamda1-lamda2  <= %20.10f\n",lamdaL1-lamda2,lamdaU1-lamda2);
	  c1=samplemean-samplemean2+Xc1[1]*sigma;
      c2=samplemean-samplemean+Xc2[1]*sigma;
      printf(" 99%% C.I. for lamda1-lamda2\n");
      lamda2=Continuous_Bernoulli_pointer(samplemean2);
      lamdaL1=Continuous_Bernoulli_pointer(samplemean+Xc1[2]*sigma);
      lamdaU1=Continuous_Bernoulli_pointer(samplemean+Xc2[2]*sigma);
      if (lamda2<0 || lamdaL1<0 || lamdaU1<0)
	   { printf(" The  99%% C.I. for lamda1-lamda2 cannot be done\n, this program wiil be terminated.");
	     goto end;}
	  printf(" %20.10f<=  lamda1-lamda2  <= %20.10f\n",lamdaL1-lamda2,lamdaU1-lamda2);
  end:
  	  system("pause");
	  free(Xc1); free(Xc2); free(filename); free(ss);
	  return 0;
}
      
      
