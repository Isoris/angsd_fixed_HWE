# ANGSD local patch notes for Hobs / hetFreq and safer site filtering

This patch set makes three focused fixes:

1. **`abcHWE` gets a new `-doHetFreq 1` mode** so you can print `hetFreq` in `.hwe.gz` **without** activating het-frequency filtering.
2. **`abcFilterSNP` no longer assumes HWE results always exist**, so SNPStat can print `NA:NA` instead of crashing or dereferencing a null pointer.
3. **`abcFilter` now bounds-checks `p->posi[s]` before indexing `fl->keeps[pos]`**, preventing out-of-bounds reads and random keep/drop behavior.

These changes are specifically aimed at safely extracting Claire Mérot-style exploratory heterozygosity (`hetFreq`) while keeping HWE/het filtering disabled by default, and at preventing hidden failures in `-sites`-based filtering.

---

## Files changed

- `abcHWE.h`
- `abcHWE.cpp`
- `abcFilter.cpp`
- `abcFilterSNP.cpp`

---

## Why each change was made

### 1) `abcHWE`: separate output from filtering
Before patching, `hetFreq` was only printed when `-minHetFreq` or `-maxHetFreq` were set, and those same flags also filtered sites. That meant there was no clean “output only” mode.

Now there is a new option:

```cpp
-doHetFreq 1
```

This prints the `hetFreq` column in `.hwe.gz` **without filtering any sites**, unless `-minHetFreq` or `-maxHetFreq` are explicitly set.

### 2) `abcFilterSNP`: guard missing HWE output
Before patching, SNPStat unconditionally dereferenced `pars->extras[9]` as `funkyHWE*`. If HWE was absent, that could crash or behave unpredictably.

Now it checks whether HWE results exist. If not, it prints:

```cpp
NA:NA
```

for `HWE_LRT:HWE_pval`, and continues the other SNPStat tests.

### 3) `abcFilter`: bounds-check `-sites` lookup
Before patching, the code directly did:

```cpp
fl->keeps[p->posi[s]]
```

with no bounds check. If `p->posi[s]` was outside `[0, target_len)`, that could read invalid memory and cause crashes or random filtering.

Now it checks:

```cpp
if(pos < 0 || pos >= L) {
  keepSites[s] = 0;
  continue;
}
```

---

## Verbatim code changes

### `abcHWE.h`

```cpp
#include "abc.h"


typedef struct {
  double *freq;
  double *freqHWE;
  double *like0;
  double *likeF;
  double *pval;
  
  double *F;

  double *hetfreq;//nspope;hetFilter
}funkyHWE;


class abcHWE:public abc{
private:
  Chisqdist *chisq;
  //  int doSnpStat;
  BGZF* outfileZ;
  kstring_t bufstr;
  void estHWE(double *x,double *loglike,int nInd);
  double HWE_like(double *x,double *loglike,int nInd);
  void HWE_EM(double *x,double *loglike,int nInd);
  //  double HWE_pval;
  int doHWE;
  double minHWEpval;
  double maxHWEpval;
  double maxHetFreq;//nspope;hetFilter
  double minHetFreq;//nspope;hetFilter
  int doHetFreq; // output hetFreq column without filtering
  int testMe;
  double tolStop;
  double differ(double *x,double *y);

public:
  abcHWE(const char *outfiles,argStruct *arguments,int inputtype);
  ~abcHWE();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);
  void print(funkyPars *pars);
  void printArg(FILE *argFile);

};
```

### `abcHWE.cpp`

```cpp
/*

  test for HWE from genotype likelihhoods


  Anders albrecht@binf.ku.dk 10 April 2016

  Authors of this file:
  Anders 

  part of ANGSD: http://www.popgen.dk/angsd 

  

*/


#include <htslib/kstring.h>
#include "abcFreq.h"
#include "abcHWE.h"
#include "aio.h"

void abcHWE::printArg(FILE *argFile){
  fprintf(argFile,"-------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-doHWE\t%d\n",doHWE);
  fprintf(argFile,"\t-minHWEpval\t%f\n",minHWEpval);
  fprintf(argFile,"\t-maxHWEpval\t%f\n",maxHWEpval);
  fprintf(argFile,"\t-maxHetFreq\t%f\n",maxHetFreq);//nspope;hetFilter
  fprintf(argFile,"\t-minHetFreq\t%f\n",minHetFreq);//nspope;hetFilter
  fprintf(argFile,"\t-doHetFreq\t%d\n",doHetFreq);
  fprintf(argFile,"\n");
}



void abcHWE::getOptions(argStruct *arguments){
  doHWE = angsd::getArg("-doHWE",doHWE,arguments);
  if(doHWE==0)
    return ;
  int GL=0;
  int doMajorMinor=0;

  minHWEpval = angsd::getArg("-minHWEpval",minHWEpval,arguments);
  maxHWEpval = angsd::getArg("-maxHWEpval",maxHWEpval,arguments);
  maxHetFreq = angsd::getArg("-maxHetFreq",maxHetFreq,arguments);//nspope;hetFilter
  minHetFreq = angsd::getArg("-minHetFreq",minHetFreq,arguments);//nspope;hetFilter
  doHetFreq = angsd::getArg("-doHetFreq",doHetFreq,arguments);
  doMajorMinor=angsd::getArg("-doMajorMinor",doMajorMinor,arguments);
  GL=angsd::getArg("-GL",GL,arguments);
  if(arguments->inputtype!=INPUT_VCF_GL && arguments->inputtype!=INPUT_GLF && arguments->inputtype!=INPUT_GLF3 && GL==0){
    fprintf(stderr,"\t-> Potential problem: You are required to choose a genotype likelihood model (-GL) \n");
    exit(0);
  } 

  if(doMajorMinor==0 && arguments->inputtype!=INPUT_VCF_GL){
      fprintf(stderr,"\t-> Potential problem: You are required to choose a major/minor estimator (-doMajorMinor)\n");
      exit(0);
  } 

  if(arguments->inputtype==INPUT_BEAGLE||arguments->inputtype==INPUT_VCF_GP){
    fprintf(stderr,"Error: you cannot estimate HWE based on posterior probabilities \n");
    exit(0);
  }

  chisq = new Chisqdist(1);//<- 1degree of freedom;

}

abcHWE::abcHWE(const char *outfiles,argStruct *arguments,int inputtype){
  chisq = NULL;
  outfileZ = NULL;

  doHWE = 0;
  maxHWEpval = -1;
  minHWEpval = -1;
  maxHetFreq = -1;//nspope;hetFilter
  minHetFreq = -1;//nspope;hetFilter
  doHetFreq = 0;
  testMe=0;
  tolStop = 0.00001;
  bufstr.s=NULL;bufstr.l=bufstr.m=0;

  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-doHWE")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }

  getOptions(arguments);
  if(doHWE==0)
    return;
  
  printArg(arguments->argumentFile);
  const char* postfix;
  postfix=".hwe.gz";
  outfileZ = aio::openFileBG(outfiles,postfix);

  //print header
  if (doHetFreq == 1 || maxHetFreq != -1 || minHetFreq != -1) //nspope;maxHetFilter
  {
    const char *str = "Chromo\tPosition\tMajor\tMinor\thweFreq\tFreq\tF\tLRT\tp-value\thetFreq\n";
    aio::bgzf_write(outfileZ,str,strlen(str));
  } else {
    const char *str = "Chromo\tPosition\tMajor\tMinor\thweFreq\tFreq\tF\tLRT\tp-value\n";
    aio::bgzf_write(outfileZ,str,strlen(str));
  }
  
}


abcHWE::~abcHWE(){

  if(outfileZ!=NULL)
    bgzf_close(outfileZ);
  delete chisq;
}


void abcHWE::clean(funkyPars *pars){
  if(doHWE==0)
    return;

  funkyHWE *hweStruct =(funkyHWE *) pars->extras[index];
  delete[] hweStruct->freq;
  delete[] hweStruct->freqHWE;
  delete[] hweStruct->hetfreq;//nspope;hetFilter
  delete[] hweStruct->F;
  delete[] hweStruct->like0;
  delete[] hweStruct->pval;
  delete[] hweStruct->likeF;
  delete hweStruct;
  
}

void abcHWE::print(funkyPars *pars){
  if(doHWE<=0)
    return;

  funkyHWE *hweStruct = (funkyHWE *) pars->extras[index];//new
  bufstr.l=0;
  for(int s=0;s<pars->numSites;s++){
    if(pars->keepSites[s]==0) 
      continue;
   
    float lrt= 2*hweStruct->like0[s]-2*hweStruct->likeF[s];
  
    if (doHetFreq==1 || maxHetFreq!=-1 || minHetFreq!=-1){//nspope;hetFilter
      ksprintf(&bufstr,"%s\t%d\t%c\t%c\t%f\t%f\t%f\t%e\t%e\t%f\n",header->target_name[pars->refId],pars->posi[s]+1,intToRef[pars->major[s]],intToRef[pars->minor[s]],hweStruct->freqHWE[s],hweStruct->freq[s],hweStruct->F[s],lrt,hweStruct->pval[s],hweStruct->hetfreq[s]);
    } else {
      ksprintf(&bufstr,"%s\t%d\t%c\t%c\t%f\t%f\t%f\t%e\t%e\n",header->target_name[pars->refId],pars->posi[s]+1,intToRef[pars->major[s]],intToRef[pars->minor[s]],hweStruct->freqHWE[s],hweStruct->freq[s],hweStruct->F[s],lrt,hweStruct->pval[s]);
    }

  }
  aio::bgzf_write(outfileZ,bufstr.s,bufstr.l);bufstr.l=0;
}


void abcHWE::run(funkyPars *pars){
 
  if(doHWE==0)
    return;

  funkyHWE *hweStruct = new funkyHWE;

  double *freq = new double[pars->numSites];
  double *freqHWE = new double[pars->numSites];
  double *hetfreq = new double[pars->numSites];//nspope;hetFilter
  double *F = new double[pars->numSites];
  double *like0 = new double[pars->numSites];
  double *likeF = new double[pars->numSites];
  double *pval = new double[pars->numSites];
   
  double **loglike3;
  loglike3=angsd::get3likesRescale(pars);


  for(int s=0;s<pars->numSites;s++){
    if(pars->keepSites[s]==0) 
      continue;
    //est under HWE
    freqHWE[s] = angsd::estFreq(loglike3[s],pars->nInd);



    //start parameter
    double x[2];
    x[1]=0.0; //F
    x[0]=freqHWE[s]; // freq

    //log like for HWE freq
    like0[s] = HWE_like(x,loglike3[s],pars->nInd);


    //start parameter for EM_F
    x[0]=freqHWE[s];
    x[1]=0.05;

    estHWE(x,loglike3[s],pars->nInd);
    freq[s]=x[0];
    F[s]=x[1];
    likeF[s] = HWE_like(x,loglike3[s],pars->nInd);

    float lrt= 2*like0[s]-2*likeF[s];
    //DRAGON lrt is sometimes nan.
    // AA: have not observed it
    pval[s]=-1;
    if(std::isnan(lrt)){
      pval[s]=lrt;
      fprintf(stdout,"nan in HWE - skipping\t %f\t%f\t%f\n",lrt,like0[s],likeF[s]);
      continue;
    }
    else if(lrt<0)
      pval[s]=1;
    else			
      pval[s]=1-chisq->cdf(lrt);

    if(maxHWEpval!=-1 && pval[s] > maxHWEpval)
      pars->keepSites[s] = 0;
    
    if(minHWEpval!=-1 && pval[s] < minHWEpval)
      pars->keepSites[s] = 0;

    //nspope;hetFilter 
    hetfreq[s] = 2.*(1.-x[1])*x[0]*(1.-x[0]);
    if( (maxHetFreq!=-1 && hetfreq[s] > maxHetFreq) || 
        (minHetFreq!=-1 && hetfreq[s] < minHetFreq) ) 
      pars->keepSites[s] = 0;
    
  }


  hweStruct->freq=freq;
  hweStruct->freqHWE=freqHWE;
  hweStruct->hetfreq=hetfreq;//nspope;hetFilter
  hweStruct->F=F;
  hweStruct->like0=like0;
  hweStruct->pval=pval;
  hweStruct->likeF=likeF;
  pars->extras[index] = hweStruct;


  for(int s=0;s<pars->numSites;s++)
    delete[] loglike3[s];
  delete[] loglike3;

}

double abcHWE::differ(double *x,double *y){
  return (fabs(x[1]-y[1]) + fabs(x[0]-y[0]));
}


void abcHWE::HWE_EM(double *x,double *loglike,int nInd){
  // freq2HWE <- function(x) {  freq <- x[2]*0.5 + x[1];   F <-	1 - x[2]/(2*freq*(1-freq)); c(freq,F);}


  

  double freq0,freq1,freq2;
  if(1){ // use F formulation instead of p0,p1,p2: compatable with the likelihood function
    double freq = x[0];
    double Fadj = freq*(1-freq)*x[1];
    freq0= pow(1-freq,2) + Fadj;
    freq1= 2*freq*(1-freq) - 2*Fadj;
    freq2= pow(freq,2) + Fadj;
  }
  else{
    freq2 = x[0];
    freq1 = x[1];
    freq0 = 1- x[0] - x[1];
  }
 
  double newFreq0=0;
  double newFreq1=0;
  double newFreq2=0;
  double norm;

  for(int i=0;i<nInd;i++){
    norm=angsd::addProtect3(log(freq0)+loglike[i*3+0],log(freq1)+loglike[i*3+1],log(freq2)+loglike[i*3+2]);
    newFreq0 += exp(log(freq0)+loglike[i*3+0] - norm);
    newFreq1 += exp(log(freq1)+loglike[i*3+1] - norm);
    newFreq2 += exp(log(freq2)+loglike[i*3+2] - norm);
  }

  newFreq0 = newFreq0/nInd;
  newFreq1 = newFreq1/nInd;
  newFreq2 = newFreq2/nInd;

  if(newFreq0 < 0.0001)
    newFreq0 = 0.0001;

  if(newFreq1 < 0.0001)
    newFreq1 = 0.0001;

  if(newFreq2 < 0.0001)
    newFreq2 = 0.0001;

  norm = newFreq0 + newFreq1 + newFreq2;
  newFreq0 = newFreq0/norm;
  newFreq1 = newFreq1/norm;
  newFreq2 = newFreq2/norm;

  
  if(newFreq0 + newFreq1 >1.00001){
    fprintf(stderr,"something is wrong i HWE\t freq %f\t F %f\n",x[0],x[1]);
    fflush(stderr);
    exit(0);
  }
 
 
  if(1){
    if(std::isnan(newFreq1*0.5 + newFreq2)){
      fprintf(stderr,"problems in abcHWE: 	x[0] %f	x[1] %f	old 1 %f 2 %f 3 %f	new 1 %f 2 %f 3 %f \n",x[0],x[1],freq0,freq1,freq2,newFreq0,newFreq1,newFreq2);
      //      exit(0);
    }
    x[0] = newFreq1*0.5 + newFreq2;
    x[1] = 1 - newFreq1/(2*x[0]*(1-x[0])) ;
  }
  else{
    x[0] = newFreq2;
    x[1] = newFreq1;
  }


 
}



double abcHWE::HWE_like(double *x,double *loglike,int nInd){
  double freq=x[0];
  double F=x[1];
  double p0=(pow(1-freq,2)+freq*(1-freq)*F);
  double p1=(2*freq*(1-freq)-2*freq*(1-freq)*F);
  double p2=(pow(freq,2)+freq*(1-freq)*F);
  double totalLogLike=0;
  for(int i=0;i<nInd;i++)
    totalLogLike+=angsd::addProtect3(log(p0)+loglike[i*3+0],log(p1) + loglike[i*3+1],log(p2) + loglike[i*3+2]);
  return -totalLogLike;
}

void abcHWE::estHWE(double *x,double *loglike,int nInd){

  double y[2];
  y[0] = x[0]+2;
  y[1] = x[1]+2;
  int iter=50;
  double dif;
  //fprintf(stderr,"%f\n",differ(x,y));
  int printer=0;
  for(int i=0;i<iter;i++){
    HWE_EM(x,loglike,nInd);
    dif = differ(x,y);
    // fprintf(stderr,"%f\n",dif);
    if(dif<tolStop)
      break;
    y[0] = x[0];
    y[1] = x[1];

  }

}
```

### `abcFilter.cpp`

```cpp
/*
  Thorfinn 31oct 2014

  refactored version of old code.
  indexing and reading of binary representation is now in prep_sites.cpp

  Code can be optimized by skipping rest of filereading for remainer of chr
osome if we have reached the last position from the keep list.

Maybe there are some memleaks. This will have to be fixed later.
*/

#include <sys/stat.h>
#include "shared.h"
#include "abc.h"
#include "analysisFunction.h"
#include "abcFilter.h"

abcFilter::~abcFilter(){
  if(fl!=NULL){
    filt_dalloc(fl);
  }
  free(fname);
}


abcFilter::abcFilter(argStruct *arguments){
    //below if shared for all analysis classes
  strict =1;
  shouldRun[index] = 1;
  header = arguments->hd;
  revMap = arguments->revMap;
  setMinIndDepth=1;
  fl = NULL;
  //his is used by this class
  keepsChr = NULL;
  curChr = -1;
  fp = NULL;
  minInd = 0;
  fname = NULL;
  doMajorMinor =0;
  capDepth = -1;
  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-sites")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }

  //get options and print them
  getOptions(arguments);
  printArg(arguments->argumentFile);

}


void abcFilter::printArg(FILE *argFile){
  fprintf(argFile,"--------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-sites\t\t%s\t(File containing sites to keep (chr pos))\n",fname);
  fprintf(argFile,"\t-sites\t\t%s\t(File containing sites to keep (chr regStart regStop))\n",fname);
  fprintf(argFile,"\t-sites\t\t%s\t(File containing sites to keep (chr pos major minor))\n",fname);
  fprintf(argFile,"\t-sites\t\t%s\t(File containing sites to keep (chr pos major minor af ac an))\n",fname);
  fprintf(argFile,"\t-minInd\t\t%d\tOnly use site if atleast minInd of samples has data\n",minInd);
  fprintf(argFile,"\t-setMinDepthInd\t%d\tOnly use site if atleast minInd of samples has this minimum depth \n",setMinIndDepth);
  fprintf(argFile,"\t-capDepth\t%d\tOnly use the first capDepth bases\n",capDepth);
  fprintf(argFile,"\t-strict\t%d\t (experimental)\n",strict);
  fprintf(argFile,"\t1) You can force major/minor by -doMajorMinor 3\n\tAnd make sure file contains 4 columns (chr tab pos tab major tab minor)\n");
}

void abcFilter::getOptions(argStruct *arguments){
  fname=angsd::getArg("-sites",fname,arguments);
  if(fname!=NULL)  
    fl = filt_read(fname);
  if(fl!=NULL)
    fprintf(stderr,"\t-> [%s] -sites is still beta, use at own risk...\n",__FILE__);


  //1=bim 2=keep
  doMajorMinor = angsd::getArg("-doMajorMinor",doMajorMinor,arguments);
  if(doMajorMinor==3 && fl!=NULL&& fl->hasExtra!=1){
    fprintf(stderr,"\t-> Must supply -sites with a file containing major and minor if -doMajorMinor 3\n");
  }
  if(doMajorMinor!=3 && fl!=NULL&& fl->hasExtra==1){
    fprintf(stderr,"\t-> Filter file contains major/minor information to use these in analysis supper '-doMajorMinor 3'\n");
  }
  capDepth = angsd::getArg("-capDepth",capDepth,arguments);
  minInd = angsd::getArg("-minInd",minInd,arguments);
  setMinIndDepth = angsd::getArg("-setMinDepthInd",setMinIndDepth,arguments);
  strict = angsd::getArg("-strict",strict,arguments);
  if(minInd >arguments->nInd){
    fprintf(stderr,"\t-> Potential problem you  filter -minInd %d but you only have %d samples?\n",minInd,arguments->nInd);
    exit(0);
  }
}


void abcFilter::run(funkyPars *p){
  //if we are comping form BAM/CRAM then we dont need to set the below
  if(p->keepSites==NULL){
    p->keepSites=new int[p->numSites];
    for(int s=0;s<p->numSites;s++)
      p->keepSites[s]=p->nInd;
  }
  if(capDepth!=-1){
    for(int s=0;s<p->numSites;s++)
      for(int i=0;i<p->nInd;i++)
	if(p->chk->nd[s][i]&&p->chk->nd[s][i]->l>capDepth)
	  p->chk->nd[s][i]->l=capDepth;
  }

  if(fl!=NULL) {
    if(fl->keeps==NULL){
      for(int s=0;s<p->numSites;s++)
	p->keepSites[s] =0;
    }else{
      for(int s=0;s<p->numSites;s++){
	//	fprintf(stderr,"p->posi:%d keeppre:%d\n",p->posi[s]+1,p->keepSites[s]);
	int pos = p->posi[s];
	int L = header->target_len[p->refId];
	if(pos < 0 || pos >= L){
	  p->keepSites[s] = 0;
	  continue;
	}
	if(strict && fl->keeps[pos]==0){
	  //	  fprintf(stderr,"never here\n");
	  p->keepSites[s] =0;
	}
	//fprintf(stderr,"p->posi:%d keeppost:%d\n",p->posi[s]+1,p->keepSites[s]);
      }
    }
  }
  //how set the keepsites according the effective sample size persite
  //if(0!=minInd){
    if(p->chk!=NULL){
      //loop over sites;
      for(int s=0;s<p->numSites;s++){
	if(p->keepSites[s]==0)
	  continue;
	int nInfo =0;
	tNode **tn = p->chk->nd[s];
	//loop over samples;
	for(int i=0;i<p->nInd;i++){
	  if(tn[i]&&tn[i]->l>=setMinIndDepth)
	    nInfo++;
	}
	p->keepSites[s] =nInfo;
	if(0!=minInd){
	  if(minInd>nInfo)
	    p->keepSites[s] =0;
	}
	//	fprintf(stderr,"p->posi:%d keeppost2:%d\n",p->posi[s]+1,p->keepSites[s]);
    }
  }
}
void abcFilter::print(funkyPars *p){
}

void abcFilter::clean(funkyPars *p){
  
}


void abcFilter::readSites(int refId) {
  //fprintf(stderr,"[%s].%s():%d -> refId:%d\n",__FILE__,__FUNCTION__,__LINE__,refId);
  if(fl==NULL)
    return;
  filt_readSites(fl,header->target_name[refId],header->target_len[refId]);
}
```

### `abcFilterSNP.cpp`

```cpp
/*
  little class that does 
1) hwe using genotype likelihoods
2) a) fisher exact
   b) GATK approach
   c) SB
   These are based on guo 2012 BMC genomics 2013 13:666
3) 2 sample wilcox rank sum for baseQ bias
   This file should be completely rewritten, much to slow and stupid. But it should give correct results

 */


#include <cmath>
#include <ctype.h>
#include "analysisFunction.h"
#include "shared.h"
#include "fet.c"
#include "chisquare.h" //<- stuff related to calculating pvalues from LRT tests
#include "abc.h"
#include "abcFilterSNP.h"
#include "abcHWE.h"
#include "abcCallGenotypes.h"
#include <htslib/kstring.h>
#include "aio.h"

/*
  wilcox manwhitney u test, whatever ,implementation is from wiki
  validated with wilcox.test (qscores of major,qscores of minor,exact=T,correct=F)
  returns Zscore, for some reason...
 */

double mann(int majD[255],int minD[255]){
  double U[255];
  double cumLess1 =0;
  for(int i=0;i<255;i++){
    U[i]=(double) minD[i]*(majD[i]/2.0+cumLess1);
    cumLess1 += majD[i];
  }
  double U1 =0;
  for(int i=0;i<255;i++) 
    U1+= U[i];
  //below is a check
  double cumLess2 =0;
  for(int i=0;i<255;i++){
    U[i]=(double) majD[i]*(minD[i]/2.0+cumLess2);
    cumLess2 += minD[i];
  }
  double U2 =0;
  for(int i=0;i<255;i++) 
    U2+= U[i];
  

  double mu=cumLess1*cumLess2/2.0;
  double sigu=sqrt((cumLess1*cumLess2*(cumLess1+cumLess2+1))/12.0);
  double Z=(std::min(U1,U2)-mu)/sigu;

 
  //  fprintf(stderr,"U1:%f U2:%f U1+U2:%f nObs:%f Z:%f p.value:%f\n",U1,U2,U1+U2,cumLess1*cumLess2,Z,2*phi(Z));

  return Z;

}



double edgebias(tNode **tn,int nInd,int maj,int min){
  //  fprintf(stderr,"phi:%f\n",phi(3));
  int majD[255];
  int minD[255];
  memset(majD,0,sizeof(int)*255);
  memset(minD,0,sizeof(int)*255);

  for(int i=0;i<nInd;i++){
    tNode *nd = tn[i];
    if(nd==NULL)
      continue;
    for(int l=0;l<nd->l;l++){
      int obB = refToInt[nd->seq[l]];

      if(obB==maj){
	majD[std::min(nd->posi[l],nd->isop[l])]++;
	//	fprintf(stdout,"maj\t%d\n",nd->qs[l]);
      }else if(obB==min){
	minD[std::min(nd->isop[l],nd->posi[l])]++;
	//fprintf(stdout,"min\t%d\n",nd->qs[l]);
      }
    }
  }

  return mann(majD,minD);
}



double mapQbias(tNode **tn,int nInd,int maj,int min){
  //  fprintf(stderr,"phi:%f\n",phi(3));
  int majD[255];
  int minD[255];
  memset(majD,0,sizeof(int)*255);
  memset(minD,0,sizeof(int)*255);

  for(int i=0;i<nInd;i++){
    tNode *nd = tn[i];
    if(nd==NULL)
      continue;
    for(int l=0;l<nd->l;l++){
      int obB = refToInt[nd->seq[l]];

      if(obB==maj){
	majD[nd->mapQ[l]]++;
	//	fprintf(stdout,"maj\t%d\n",nd->qs[l]);
      }else if(obB==min){
	minD[nd->mapQ[l]]++;
	//fprintf(stdout,"min\t%d\n",nd->qs[l]);
      }
    }
  }

  return mann(majD,minD);
}



double baseQbias(tNode **tn,int nInd,int maj,int min){
  //  fprintf(stderr,"phi:%f\n",phi(3));
  int majD[255];
  int minD[255];
  memset(majD,0,sizeof(int)*255);
  memset(minD,0,sizeof(int)*255);

  for(int i=0;i<nInd;i++){
    tNode *nd = tn[i];
    if(nd==NULL)
      continue;
    for(int l=0;l<nd->l;l++){
      int obB = refToInt[nd->seq[l]];

      if(obB==maj){
	majD[nd->qs[l]]++;
	//	fprintf(stdout,"maj\t%d\n",nd->qs[l]);
      }else if(obB==min){
	minD[nd->qs[l]]++;
	//fprintf(stdout,"min\t%d\n",nd->qs[l]);
      }
    }
  }

  return mann(majD,minD);
}

Chisqdist chi(1);

//guo 2012 mutat res 2012, 744(2):154-160
double sb1(int cnts[4]){
  double a=cnts[0];double b=cnts[1];double c=cnts[2];double d=cnts[3];
  double top=b/(a+b)-d/(c+d);
  double bot=(b+d)/(a+b+c+d);
  return top/bot;
}

//the gatk way
double sb2(int cnts[4]){
  double a=cnts[0];double b=cnts[1];double c=cnts[2];double d=cnts[3];
  double en=(b/(a+b))*(c/(c+d));
  double to=(a+c)/(a+b+c+d);
  double tre=(d/(c+d))*(a/(a+b));
  return std::max(en/to,tre/to);
}

//strandbias using fisher

double sb3(int cnts[4]){

  double left,right,twotail,prob;
  kt_fisher_exact(cnts[0], cnts[1], cnts[2], cnts[3], &left, &right, &twotail);
  return twotail;
}


void abcFilterSNP::printArg(FILE *argFile){
   fprintf(argFile,"-----BETA---------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-doSnpStat %d\n",doSnpStat);
  fprintf(argFile,"\t-edge_pval %f\n",edge_pval);
  fprintf(argFile,"\t-mapQ_pval %f\n",mapQ_pval);
  fprintf(argFile,"\t-sb_pval %f\n",sb_pval);
  fprintf(argFile,"\t-hwe_pval %f\n",hwe_pval);
  fprintf(argFile,"\t-qscore_pval %f\n",qscore_pval);
  fprintf(argFile,"\t-hetbias_pval %f\n",hetbias_pval);

}
void abcFilterSNP::run(funkyPars *pars){
  if(!doSnpStat)
    return;
  chunkyT *chk = pars->chk;
  
  if(doSnpStat==1){
    kstring_t *bufstr = new kstring_t;
    bufstr->s=NULL;bufstr->l=bufstr->m=0;
    //loop over sites;
    kstring_t persite;
    persite.s=NULL;persite.l=persite.m=0;
    //pull 
   
    for(int s=0;s<pars->numSites;s++) {
      if(pars->keepSites[s]==0)
	continue;
      //loop over samples
      int cnts[4]={0,0,0,0};
      for(int i=0;i<pars->nInd;i++){
	tNode *nd = chk->nd[s][i];
	if(nd==NULL)
	  continue;
	for(int l=0;l<nd->l;l++){
	  int obB = refToInt[nd->seq[l]];
	  //	    fprintf(stderr,"%c ",nd.seq[l]);
	  int strand = (isupper(nd->seq[l])==0)<<1;
	  //  fprintf(stderr,"strand:%d\n",strand);
	  if(obB==4)
	    continue;
	  if((obB!=pars->major[s] && obB!=pars->minor[s]) )
	    continue;
	  if(obB!=pars->major[s])
	    strand +=1;
	  //fprintf(stderr,"strand=%d\n",strand);
	  cnts[strand]++;
	}
      }
      ksprintf(&persite,"%s\t%d\t%d %d %d %d\t",header->target_name[pars->refId],pars->posi[s]+1, cnts[0],cnts[1],cnts[2],cnts[3]);
      ksprintf(&persite,"%f:%f:%f\t",sb1(cnts),sb2(cnts),sb3(cnts));
      //      fprintf(stderr,"sb4:%f\n",sb3(cnts));
      if(sb_pval!=-1 && sb3(cnts)<sb_pval)
	pars->keepSites[s] = 0;

      funkyHWE *hweStruct = (funkyHWE *) pars->extras[9];//THIS IS VERY NASTY! the ordering within general.cpp is now important
      if(hweStruct){
	double lrt = 2*hweStruct->like0[s]-2*hweStruct->likeF[s];
	double pval;
	if(std::isnan(lrt))
	  pval=lrt;
	else if(lrt<0)
	  pval =1;
	else
	  pval =1- chi.cdf(lrt);
	ksprintf(&persite,"%f:%e\t",lrt,pval);
	if(hwe_pval!=-1 && !std::isnan(pval) && pval<hwe_pval)
	  pars->keepSites[s] = 0;
      } else {
	ksprintf(&persite,"NA:NA\t");
      }

      double Z = baseQbias(chk->nd[s],pars->nInd,refToInt[pars->major[s]],refToInt[pars->minor[s]]);
      ksprintf(&persite,"%f:%e\t",Z,2*phi(Z));
      if(qscore_pval!=-1 && 2*phi(Z)<qscore_pval)
	pars->keepSites[s] = 0;
    

      Z = mapQbias(chk->nd[s],pars->nInd,refToInt[pars->major[s]],refToInt[pars->minor[s]]);
      ksprintf(&persite,"%f:%e\t",Z,2*phi(Z));

      if(mapQ_pval!=-1 && 2*phi(Z)<mapQ_pval)
	pars->keepSites[s] = 0;

      Z = edgebias(chk->nd[s],pars->nInd,refToInt[pars->major[s]],refToInt[pars->minor[s]]);
      ksprintf(&persite,"%f:%e",Z,2*phi(Z));
      if(edge_pval!=-1 && 2*phi(Z)<edge_pval)
	pars->keepSites[s] = 0;

      genoCalls *gcw =(genoCalls *) pars->extras[11];
      int **gc=NULL;
      if(gcw)
	gc = gcw->dat;
      if(gc){
	cnts[0]=cnts[1]=cnts[2]=cnts[3]=0;

	int nsampleswithdata =0;
	for(int i=0;i<pars->nInd;i++){
	  if(gc[s][i]!=1)
	    continue;
	  tNode *nd = chk->nd[s][i];
	  if(nd==NULL)
	    continue;
	  nsampleswithdata++;
	  for(int l=0;l<nd->l;l++){
	    int obB = refToInt[nd->seq[l]];
	    //	    fprintf(stderr,"%c ",nd.seq[l]);
	    int strand = (isupper(nd->seq[l])==0)<<1;
	    //  fprintf(stderr,"strand:%d\n",strand);
	    if(obB==4)
	      continue;
	    if((obB!=pars->major[s] && obB!=pars->minor[s]) )
	      continue;
	    if(obB!=pars->major[s])
	      strand +=1;
	    //fprintf(stderr,"strand=%d\n",strand);
	    cnts[strand]++;
	  }
	}

	double tsum= cnts[0]+cnts[1]+cnts[2]+cnts[3];
	int n = tsum;
	ksprintf(&persite,"\t%d %d %d %d %d\t",cnts[0],cnts[1],cnts[2],cnts[3],n);

	if(tsum>0){
	  double fA=(cnts[0]+cnts[2])/tsum;
	  double fa=(cnts[1]+cnts[3])/tsum;
	  double lrt = 2*tsum*(fA-0.5)*(fA-0.5) + 2*tsum*(fa-0.5)*(fa-0.5);
	  double pval;
	  if(std::isnan(lrt))
	    pval=lrt;
	  else if(lrt<0)
	    pval =1;
	  else
	    pval =1- chi.cdf(lrt);
	  ksprintf(&persite,"%f:%e\t",lrt,pval);
	  if(hetbias_pval!=-1 && !std::isnan(pval) && pval<hetbias_pval)
	    pars->keepSites[s]=0;
	} else {
	  ksprintf(&persite,"nan:nan\t");
	}
      }

      
      if(pars->keepSites[s]!=0&&persite.l>0){
	ksprintf(&persite,"\n");
	ksprintf(bufstr,"%s",persite.s);
      }
      persite.l=0;
   }
    free(persite.s);
    pars->extras[index] = bufstr;
  }
}

void abcFilterSNP::clean(funkyPars *fp){
  if(!doSnpStat)
    return;

  
}

void abcFilterSNP::print(funkyPars *pars){
  if(!doSnpStat)
    return;
  kstring_t *bufstr =(kstring_t*) pars->extras[index];
  aio::bgzf_write(outfileZ,bufstr->s,bufstr->l);bufstr->l=0;
  free(bufstr->s);
  delete bufstr;
}


void abcFilterSNP::getOptions(argStruct *arguments){
  //default


  //from command line
  doSnpStat=angsd::getArg("-doSnpStat",doSnpStat,arguments);


  if(doSnpStat==0)
    return;
  int domajorminor=0;
  domajorminor = angsd::getArg("-domajorminor",domajorminor,arguments);
  if(domajorminor==0){
    fprintf(stderr,"\t-> Must supply -doMajorMinor for running dosnpstat (needs to look a distributions of major and minor alleles)\n");
    exit(0);
  }
  //from command line

  edge_pval=angsd::getArg("-edge_pval",edge_pval,arguments);      
  mapQ_pval=angsd::getArg("-mapQ_pval",mapQ_pval,arguments);      
  sb_pval=angsd::getArg("-sb_pval",sb_pval,arguments);    
  hwe_pval=angsd::getArg("-hwe_pval",hwe_pval,arguments);    
  qscore_pval=angsd::getArg("-qscore_pval",qscore_pval,arguments);
  hetbias_pval=angsd::getArg("-hetbias_pval",hetbias_pval,arguments);    
  int doHWE=0;
  doHWE=angsd::getArg("-doHWE",doHWE,arguments);    
  if(doHWE==0){
    fprintf(stderr,"must use -doHWE 1, to test for HWE\n");
    exit(0);
  }  
}


abcFilterSNP::abcFilterSNP(const char *outfiles,argStruct *arguments,int inputtype){
  doSnpStat=0;
  outfileZ = NULL;
  edge_pval = -1;
  mapQ_pval = -1;
  sb_pval = -1;
  hwe_pval = -1;
  qscore_pval = -1;
  hetbias_pval = -1;
  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-doSnpStat")||!strcasecmp(arguments->argv[1],"-doPost")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }

  getOptions(arguments);
  
  if(doSnpStat==0){
    shouldRun[index] =0;
    return;
  }
  printArg(arguments->argumentFile);  
  int dogeno=0;
  dogeno=angsd::getArg("-dogeno",dogeno,arguments);
  if(doSnpStat){
    //    fprintf(stderr,"running doSnpStat=%d\n",doSnpStat);
    const char *postfix=".snpStat.gz";
    outfileZ = aio::openFileBG(outfiles,postfix);
    kstring_t bufstr;bufstr.s=NULL;bufstr.l=bufstr.m=0;
    ksprintf(&bufstr,"Chromo\tPosition\t+Major +Minor -Major -Minor\tSB1:SB2:SB3\tHWE_LRT:HWE_pval\tbaseQ_Z:baseQ_pval\tmapQ_Z:mapQ_pval\tedge_z:edge_pval");
    if(dogeno){
      ksprintf(&bufstr,"\t+MajorHet +MinorHet -MajorHet -MinorHet nHet\thetStat:hetStat_pval");

    }
    ksprintf(&bufstr,"\n");
    aio::bgzf_write(outfileZ,bufstr.s,bufstr.l);bufstr.l=0;
    free(bufstr.s);
  }

}

abcFilterSNP::~abcFilterSNP(){
  if(outfileZ!=NULL)
    bgzf_close(outfileZ);

}
```

---

## Recommended CLI after patching

For exploratory Claire Mérot-style heterozygosity output, use:

```bash
angsd \
  -b bamlist.txt \
  -ref ref.fa \
  -r chr1 \
  -GL 1 \
  -doMajorMinor 1 \
  -doMaf 1 \
  -doHWE 1 \
  -doHetFreq 1 \
  -out test_chr1
```

Do **not** set `-minHetFreq`, `-maxHetFreq`, `-minHWEpval`, or `-maxHWEpval` unless you intentionally want filtering.

---

## Rebuild

Typical rebuild:

```bash
make clean
make -j
```

---

## Quick validation

After running a test chromosome:

```bash
zcat test_chr1.hwe.gz | head
```

You want to see:
- a `hetFreq` column
- normal site counts
- no massive silent collapse of sites

---

## Notes

- `theta/pi` tracks remain separate from `hetFreq/Hobs`.
- `-doHetFreq 1` is an **output** switch only.
- `-minHetFreq` / `-maxHetFreq` remain **filter** switches.
- `abcFilterSNP` now tolerates missing HWE state safely.
