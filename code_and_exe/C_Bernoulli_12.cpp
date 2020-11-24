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
// computing the ANOVA test statistic sampling distribution and p value
double Continuous_Bernoulli_ANOVA(int k,int nok,double lamda,double Ftest)
{
	int n=1000000,i,ii,j,*samplesize,nT;
	double *samplemean,**sampledata,grandmean,SST,SSTR,SSE,*getdata,pvalue;
	
	  samplesize=(int *)malloc(k*sizeof(int));// sample size vector
      for(i=0;i<k;++i) samplesize[i]=nok; // the sample size is nok
	  getdata=(double *)malloc(n*sizeof(double)); // test statistic database
	  nT=0;
      for(i=0;i<k;++i) nT=nT+samplesize[i]; // total sample number
      samplemean=(double *)malloc(k*sizeof(double));  // sample mean vector
      sampledata=(double **)malloc(k*sizeof(double *)); // sample data matrix
      for(i=0;i<k;++i) 
        sampledata[i]=(double *)malloc(samplesize[i]*sizeof(double));
      for(i=0;i<n;++i)
      {
      grandmean=0; 
      for(ii=0;ii<k;++ii)
      {
	    samplemean[ii]=0;
        for(j=0;j<samplesize[ii];++j)
          {
          	sampledata[ii][j]=simulated_continue_bernoulli(lamda);
          	samplemean[ii]=samplemean[ii]+sampledata[ii][j];
          	grandmean=grandmean+sampledata[ii][j];
		  }
		samplemean[ii]=samplemean[ii]/(double)samplesize[ii];
       }
       grandmean=grandmean/nT;
       SSE=0;SST=0;
		for(ii=0;ii<k;++ii)
           for(j=0;j<samplesize[ii];++j)
           {
           	SSE=SSE+(sampledata[ii][j]-samplemean[ii])*(sampledata[ii][j]-samplemean[ii]);
           	SST=SST+(sampledata[ii][j]-grandmean)*(sampledata[ii][j]-grandmean);
		   }
		SSTR=SST-SSE;
		getdata[i]=((double)(nT-k)/(double)(k-1))*SSTR/SSE; // test statistic simulated data
	   }  
	   divide_Conquer_X_step(getdata,n); // sorting
	   for(i=0;i<n;++i) if (getdata[i]>=Ftest) break;
	   pvalue=1-(double)(i+1)/(double)n;  // p value
	   if (pvalue<0) pvalue=0;
	   if (pvalue>1) pvalue=1;
       free(getdata);
       free(samplesize);
       free(samplemean);
       for(i=0;i<k;++i) free(sampledata[i]);
       free(sampledata);
    return(pvalue);
}

int main() {  
     double lamda,*lamdaX,*samplemean,**sampledata,Ftest,R2;
     double SST,SSTR,SSE,grandmean,pvalue;
	 int i,ii,*samplesize,ret,j=0,k,nT,nok;
	 char *ss;
      
      ss=(char *)malloc(1024*sizeof(char));
      
      again1:
      printf(" treatment number, k(2~20)>");
      scanf("%s",ss);
      k=atoi(ss);
      if (k<=1 || k>20) goto again1;
      lamdaX=(double *)malloc(k*sizeof(double));
      for(i=0;i<k;++i)
      {
      	again2:
          printf(" lamda %d>",i+1);
          scanf("%s",ss);
          lamdaX[i]=atof(ss);
      if (lamdaX[i]<=0 || lamdaX[i]>=1) goto again2;
	  }
    
      again3:
      printf(" sample size of each treatment>");
      scanf("%s",ss);
      nok=atoi(ss);
      if (nok<=1) goto again3;
      samplesize=(int *)malloc(k*sizeof(int));// sample size vector
      
      for(i=0;i<k;++i) samplesize[i]=nok; // the sample size is nok
      nT=0;
      for(i=0;i<k;++i) nT=nT+samplesize[i]; // the total sample size
      samplemean=(double *)malloc(k*sizeof(double));// sample mean vector
      sampledata=(double **)malloc(k*sizeof(double *)); // sample data matrix
      for(i=0;i<k;++i) 
        sampledata[i]=(double *)malloc(samplesize[i]*sizeof(double));
      grandmean=0;
      // getting random sample from the simulator
      for(ii=0;ii<k;++ii)
      {
	    samplemean[ii]=0;
        for(j=0;j<samplesize[ii];++j)
          {
          	sampledata[ii][j]=simulated_continue_bernoulli(lamdaX[ii]);
          	samplemean[ii]=samplemean[ii]+sampledata[ii][j];
          	grandmean=grandmean+sampledata[ii][j];
		  }
		samplemean[ii]=samplemean[ii]/(double)samplesize[ii]; // each treatment sample mean
       }
       // computing the grand mean
       grandmean=grandmean/nT;
       // computing SST, SSTR, SSE
       SSE=0;SST=0;
		for(ii=0;ii<k;++ii)
           for(j=0;j<samplesize[ii];++j)
           {
           	SSE=SSE+(sampledata[ii][j]-samplemean[ii])*(sampledata[ii][j]-samplemean[ii]);
           	SST=SST+(sampledata[ii][j]-grandmean)*(sampledata[ii][j]-grandmean);
		   }
		SSTR=SST-SSE;
		Ftest=((double)(nT-k)/(double)(k-1))*SSTR/SSE;
		R2=SSTR/SST;
		for(i=0;i<k;++i) printf("X%d~Continuous Bernoulli(lamda%d=%f),\n",i+1,i+1,lamdaX[i]);
		printf(" The treatment number=%d,\n",k);
		printf(" The sample size of each treatment=%d,\n",nok);
		printf(" The sample mean of each treatment,\n");
		for(i=0;i<k;++i) printf(" Treatment %d sample mean=%20.10f,\n",i+1,samplemean[i]);
		lamda=Continuous_Bernoulli_pointer(grandmean);
		if (lamda<0) { printf(" ANOVA cannot be finished, when lamda estimate is error,\n this program will be terminated."); goto end;}
		printf(" grand mean=%20.10f, the grand lamda estimated=%20.10f\n",grandmean,lamda);
		printf(" The ANOVA,\n");
		printf(" Treatment:SSTR=%20.10f, df=%d, MSTR=%20.10f,\n",SSTR,k-1,SSTR/(double)(k-1));
		printf(" Error:    SSE =%20.10f, df=%d, MSE=%20.10f,\n",SSE,nT-k,SSE/(double)(nT-k));
		printf(" Total:    SST =%20.10f, df=%d,\n",SST,nT-1);
		printf(" Please wait a moment to get the test statistic sampling distribution,..\n");
		pvalue=Continuous_Bernoulli_ANOVA(k,nok,lamda,Ftest);
		printf(" H0:lamda1=lamda2=......lamda%d,\n",k);
		printf(" R2=%20.10f, F test=MSTR/MSE=%20.10f\n",R2,Ftest);
		printf(" p value=%f,\n",pvalue);
	    system("pause");
	    free(samplesize);
        for(i=0;i<k;++i) free(sampledata[i]);
        free(sampledata);
  end:
	  free(ss); 
	  free(samplesize);
      for(i=0;i<k;++i) free(sampledata[i]);
      free(sampledata);
	  return 0;
}
      
      
