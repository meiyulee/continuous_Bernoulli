#include <math.h>
#include <stdio.h>
#include <stdlib.h> 
#include <malloc.h>
#include <time.h> 
#include <io.h> 
void frequecy_table(int n,double *getvalue1)
{
   int i,ni,classno1,count;
   double len1,tmin1,tmax1,width1,*pdfvalue1,*upper1,*lower1,cuumlative_pr=0,class_midpoint;
   FILE *file;
   
    tmin1=getvalue1[0]; tmax1=getvalue1[0];
    for(ni=1;ni<n;++ni)
     {  if (getvalue1[ni]<tmin1) tmin1=getvalue1[ni];
        if (getvalue1[ni]>tmax1) tmax1=getvalue1[ni];
                 }
     len1=tmax1-tmin1;
     classno1=200;
     width1=len1/(double)classno1;
     pdfvalue1=(double *)malloc(classno1*sizeof(double));
     upper1=(double *)malloc(classno1*sizeof(double));
     lower1=(double *)malloc(classno1*sizeof(double));
     for(ni=0;ni<classno1;++ni)   pdfvalue1[ni]=0;   
     upper1[0]=tmin1+width1;
     for(ni=1;ni<classno1;++ni)
     {  lower1[ni]=upper1[ni-1]; upper1[ni]=lower1[ni]+width1; }
     upper1[classno1-1]=tmax1+1;
     lower1[0]=tmin1-1;
     for(ni=0;ni<n;++ni)
     { 
        count=(long)((getvalue1[ni]-tmin1)/width1);
        if (count>=0 && count<classno1) pdfvalue1[count]+=1; 
		 }
    for(ni=0;ni<classno1;++ni)  pdfvalue1[ni]/=(double)(n); 
    lower1[0]=tmin1;upper1[classno1-1]=tmax1;
    file=fopen("C:\\C_Bernoulli\\tep\\frequency_table.txt","w"); 
    cuumlative_pr=0;
    for(ni=0;ni<classno1;++ni)
    {
     class_midpoint=(lower1[ni]+upper1[ni])/2;
     cuumlative_pr=cuumlative_pr+pdfvalue1[ni];
     fprintf(file,"%20.10f          %20.10f       %20.10f         %20.10f\n",class_midpoint,pdfvalue1[ni]/width1,class_midpoint,cuumlative_pr); 
     }
   fclose(file);
   free(pdfvalue1); free(upper1); free(lower1);
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
    
    // checking the lamda is 0.5 or not 
    if (lamda>=0.499 && lamda<=0.501)
    {
    sprintf(ss,"%12.10f",lamda);
    kk=strcompare_spc(ss,(char *)"0.5000000000",12);
    	}
    // kk is 0 when lamda=0.5, X~Uniform(0,1) if lamda=0.5
    if (kk!=0) 
    {
     // the inverse distribution function is the simulated data
         tep1=simulated_uniform(0,1); 
         tep1=(2*lamda-1)*tep1-lamda+1;
         tep1=tep1/(1-lamda);
         tep=log(tep1)/log(lamda/(1-lamda));
		}
	else
    tep=simulated_uniform(0,1);  
   return(tep);
}
// computing the mean,variance,skewed coefficient and kurtosis coefficient
void moment_lamda(int n,double *datastore)
{
	int i;
	double EX,VarX,skew,kurtosis;
	double sumX,sumX2,sumX3,sumX4,sigmaX,tep;
	
	sumX=0; sumX2=0; 
    for(i=0;i<n;++i)
       	sumX=sumX+datastore[i];
	sumX=sumX/(double)n;
	for(i=0;i<n;++i)
      	sumX2=sumX2+(datastore[i]-sumX)*(datastore[i]-sumX);
	sumX2=sumX2/(double)n;
	if (sumX2>0) 
	{
		sigmaX=exp(0.5*log(sumX2)); 
		sumX3=0;
        sumX4=0;
		for(i=0;i<n;++i)
		{
		  tep=(datastore[i]-sumX)/sigmaX;
		  sumX3=sumX3+tep*tep*tep;
		 sumX4=sumX4+tep*tep*tep*tep;
			}
		sumX3=sumX3/(double)n;
		sumX4=sumX4/(double)n;
		EX=sumX;VarX=sumX2;skew=sumX3;kurtosis=sumX4;
		printf(" E(lamda hat)=%20.10f,\n Var(lamda hat)=%20.10f,\n",EX,VarX);
        printf(" skewed coefficient(lamda hat)=%20.10f,\n kurtosis coefficient(lamda hat)=%20.10f\n",skew,kurtosis);
      }
      else
       printf(" E(lamda hat)=%20.10f\n",EX);
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
int main()
{  
    int i,j,n=100000000,samplesize;
    double *datastore,lamda,tepdata,lamdahat;
    char *ss;
    
    ss=(char *)malloc(1024*sizeof(char));
    // input the parameter lamda
    loop1:
    printf(" The Continuous Bernoulli, lamda=");
	scanf("%s",ss);
	lamda=atof(ss);
	if (lamda<=0 || lamda>=1 )
	{ 
	  printf(" The Continuous Bernoulli(lamda=%f), it is error, 0<lamda<1,\n please input again.\n",lamda);
	  goto loop1;
	  }  
	loop2:
    printf(" The sample size=");
	scanf("%s",ss);
	samplesize=atoi(ss);
	if (samplesize<=0 )
	{ 
	  printf(" The sample size, it is error, sample size must greater than 0,\n please input again.\n",samplesize);
	  goto loop2;
	  }  
    datastore=(double *)malloc(n*sizeof(double));
    printf(" please wait a moment to collect the simulated data,....\n");
    // getting the simulated data and the amount=100,000,000, 
	// the moment of this database is closed to the distribution
	printf(" finished %10.8f%%\n",(double)1/(double)n*100);
	for(i=0;i<n;++i)
	{
		if ((i+1)%10000000==0) printf(" finished %10.8f%%\n",(double)(i+1)/(double)n*100);
		tepdata=0;
		for(j=0;j<samplesize;++j)
		 tepdata=tepdata+simulated_continue_bernoulli(lamda); // simulated data summation
		tepdata=tepdata/(double)samplesize; // sample mean
		lamdahat=Continuous_Bernoulli_pointer(tepdata);
		if (lamdahat<0) { --i; continue; }
		datastore[i]=lamdahat;
	}
	printf(" X1,X2,...,X%d iid~Continuous Bernoulli(lamda=%f),\n",samplesize,lamda);
	printf(" X bar=(X1+X2+....+X%d)/%d=the sample mean,\n",samplesize,samplesize);
	printf(" lamda hat is converted from X bar.\n");
	moment_lamda(n,datastore);
	printf(" the lamda hat frequency table is c:\\C_Bernoulli\\tep\\frequency_table.txt\n");
	i=chdir("c:\\C_Bernoulli"); // check this subdirectory is existed
    if (i!=0) mkdir("c:\\C_Bernoulli"); // make this subdirectory
    i=chdir("c:\\C_Bernoulli\\tep"); // check this subdirectory is existed
    if (i!=0) mkdir("c:\\C_Bernoulli\\tep"); // make this subdirectory
	frequecy_table(n,datastore);
	free(ss); free(datastore);
	system("pause");
	return 0;
}
      
      
