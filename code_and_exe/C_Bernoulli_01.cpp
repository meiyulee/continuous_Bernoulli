#include <math.h>
#include <stdio.h>
#include <stdlib.h> 
#include <malloc.h>
#include <time.h> 
#include <io.h> 
// Getting the time 
int gettimexx(void)
{
 time_t  ti;
 struct tm *timeinfo;
 int i=0;  
 
 ti=time(NULL);/*  allocate  the buffer*/
 timeinfo=localtime(&ti);  /*  catch the time*/
 i=timeinfo->tm_hour*3600+timeinfo->tm_min*60+timeinfo->tm_sec;
 return(i);
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
int main()
{  
    int loopnumber,i,number;
    double simulated_data,lamda;
    char *ss;
    FILE *file;
    
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
	// input the simulated data number
    loop2:
    printf(" The simulated data number=");
	scanf("%s",ss);
	number=atoi(ss);
	if (number<=0)
	{ 
	  printf(" The simulated data number=%d, it is error,\n please input again.\n",number);
	  goto loop2;
	  }  
	// the return time is the number of calling rand(),
	// rand() is always the same number series. This is changed the rule.
	loopnumber=gettimexx();
    for(i=0;i<loopnumber;++i) rand();
    // The simulated data 
    i=chdir("c:\\C_Bernoulli"); // check this subdirectory is existed
    if (i!=0) mkdir("c:\\C_Bernoulli"); // make this subdirectory
    i=chdir("c:\\C_Bernoulli\\tep"); // check this subdirectory is existed
    if (i!=0) mkdir("c:\\C_Bernoulli\\tep"); // make this subdirectory
    
    file=fopen("c:\\C_Bernoulli\\tep\\simulated_data.txt","w");
	for(i=0;i<number;++i)
	{
	  simulated_data=simulated_continue_bernoulli(lamda);
	  printf("%20.10f\n",simulated_data);
	  fprintf(file,"%20.10f\n",simulated_data);
	}
	fclose(file);
	printf(" the simulated data filename is c:\\C_Bernoulli\\tep\\simulated_data.txt,\n");
	free(ss);
	system("pause");
	return 0;
}
      
      
