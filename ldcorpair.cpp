#include <string.h>
#include <math.h>
#include <cstdio>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>

#define PACK_DENSITY 4
#define MASK0   3 // 3 << 2 * 0
#define MASK1  12 // 3 << 2 * 1
#define MASK2  48 // 3 << 2 * 2
#define MASK3 192 // 3 << 2 * 3

using namespace std;

void decode_plink(char *output, const char *input, const int lengthInput){
  int i, k;
  char tmp, geno;
  int a1, a2;
  
  for(i=0;i<lengthInput;++i){
    tmp = input[i];
    k   = PACK_DENSITY * i;
    geno      = (tmp & MASK0);
    a1        = !(geno & 1);
    a2        = !(geno >> 1);
    output[k] = (geno == 1) ? 3 : a1 + a2;
    k++;
    geno      = (tmp & MASK1) >> 2; 
    a1        = !(geno & 1);
    a2        = !(geno >> 1);
    output[k] = (geno == 1) ? 3 : a1 + a2;
    k++;
    geno      = (tmp & MASK2) >> 4; 
    a1        = !(geno & 1);
    a2        = !(geno >> 1);
    output[k] = (geno == 1) ? 3 : a1 + a2;
    k++;
    geno      = (tmp & MASK3) >> 6; 
    a1        = !(geno & 1);
    a2        = !(geno >> 1);
    output[k] = (geno == 1) ? 3 : a1 + a2;
  }
}

// Main funtion used on cov_files.txt with a sequential access
int main(int argc, char *argv[]){
  // Input arguments
  string prefix       = "none";
  string outfix       = "none";
  string snplistfile  = "none";
  float window_kb     =   100.f;
   
  // Indices
  string sw;
  int i,j,k;
  
  if(argc==1){
    cerr<<"\tArguments must be specified. Type --help for more details."<<endl;
    exit(1);
  }
  
  // Read arguments
  sw = argv[1];
  if (sw == "--help"){
    cerr<<"\t--bfile      : Binary PLINK format for genotypes."<<endl;
    cerr<<"\t--snplist    : Specify the list of SNPs containing the target SNPs."<<endl;
    cerr<<"\t--window-kb  : Specify the window size to calculate LD correlations. Default is 100 kb."<<endl;
    cerr<<"\t--out        : A prefix for the output file [prefix].ldcor.pairs."<<endl;
    cerr<<"\t[Note] Missing values are imputed to the major allele."<<endl;
    exit(1);
  }else{
    if (argc == 1) {
      cerr<<"\tArguments must be specified. Type --help for more details."<<endl;
      exit(1);
    }
  }
  
  for(i = 1; i<argc;i++){
    sw = argv[i];
    if (sw == "--bfile"){
      prefix = argv[i + 1];
    }
    if (sw == "--snplist"){
      snplistfile = argv[i + 1];
    }
    if (sw == "--window-kb"){
      window_kb = atof(argv[i + 1]);
    }
    if (sw == "--out"){
      outfix = argv[i + 1];
    }
  }
  
  if(prefix=="none"){
    cerr<<"\tA prefix must be specified for files [prefix].bed, [prefix].bim and [prefix].fam. Type [./ldcorpair --help.]"<<endl;
    exit(1);
  }
  if(snplistfile=="none"){
    cerr<<"\tAn input list of SNPs must be specified. Type [./ldcorpair --help.]"<<endl;
    exit(1);
  }
  if(outfix=="none"){
    cerr<<"\tAn output file prefix name must be specified. Type [./ldcorpair --help.]"<<endl;
    exit(1);
  }

  string bedfile = prefix+".bed";
  string bimfile = prefix+".bim";
  string famfile = prefix+".fam";

  // Few tools
  string line = "";
  string tok  = ""; 
  ifstream tmpStream;
  
  // Get number of SNPs
  int M = -1;
  tmpStream.open(bimfile.c_str());
  while(tmpStream){
    getline(tmpStream,line);
    M++;
  }
  tmpStream.close();
  cout<<"# Found "<<M<<" SNPs."<<endl;

  // Get sample size
  int N = -1; 
  tmpStream.open(famfile.c_str());
  while(tmpStream){
    getline(tmpStream,line);
    N++;
  }
  tmpStream.close();  
  cout<<"# Found "<<N<<" samples."<<endl;
  
  int nSNPlist = -1;
  tmpStream.open(snplistfile.c_str());
  while(tmpStream){
    getline(tmpStream,line);
    nSNPlist++;
  }
  tmpStream.close();
  
  cout<<"# Found "<<nSNPlist<<" SNPs in file named: "<<snplistfile<<".\n";
  string *snplist = new string[nSNPlist];
  tmpStream.open(snplistfile.c_str());
  for(j=0;j<nSNPlist;j++){
    tmpStream >> snplist[j];
  }
  tmpStream.close();
  
  string *CHROM  = new string[M];
  string *SNPS   = new string[M];
  int *snpsFound = new int[M];
  float  *POSKB  = new float[M];
  string *A1     = new string[M];
  string *A2     = new string[M];

  int nSNPfound = 0;
  tmpStream.open(bimfile.c_str());
  for(j=0;j<M;j++){
    tmpStream >> CHROM[j];
    tmpStream >> SNPS[j];
    tmpStream >> tok;
    tmpStream >> tok; POSKB[j] = atof(tok.c_str()) / 1000.f;
    tmpStream >> A1[j];
    tmpStream >> A2[j];
    snpsFound[j] = -1;
    for(k=0;k<nSNPlist;k++){
      if(SNPS[j]==snplist[k]){
        snpsFound[j] = k;
        nSNPfound++;
      }
    }
  }
  tmpStream.close();
  cout<<"# Found "<<nSNPfound<<" SNPs in common.\n";
  
  // Re-order indexes -- very important!!!
  string *chrom  = new string[nSNPfound];
  string *snps   = new string[nSNPfound];
  float  *poskb  = new float[nSNPfound];
  string *a1     = new string[nSNPfound];
  string *a2     = new string[nSNPfound];
  k = -1;
  for(j=0;j<M;j++){
    if(snpsFound[j] != -1){
      k++;
      snpsFound[j] = k;
      chrom[k] = CHROM[j];
      snps[k]  = SNPS[j];
      poskb[k] = POSKB[j];
      a1[k]    = A1[j];
      a2[k]    = A2[j];
    }
  }
  
  int N_times_nSNPfound = N * nSNPfound;
  double* xTarget = new double[N_times_nSNPfound];

  int numBytes   = (int)ceil((double)N / PACK_DENSITY);
  char* packed   = new char[numBytes];
  char* unpacked = new char[numBytes * PACK_DENSITY];

  // PASS 1
  cout<<"# Pass 1: loading target SNP information...";
  ifstream influx;
  influx.open(bedfile.c_str(), std::ios::in | std::ios::binary);
  if(!influx){
    cerr << "[readGenotypes] Error reading file "<<bedfile<<endl;
    exit(1);
  }
  influx.seekg(0, ifstream::end);
  influx.seekg(3, ifstream::beg);
  
  int nEff;
  for(j=0;j<M;j++){
    influx.read((char*)packed, sizeof(char) * numBytes);
    if(snpsFound[j]>=0){
      k = snpsFound[j];
      decode_plink(unpacked, packed, numBytes);
      for(i=0;i<N;i++){
        xTarget[i+k*N] = (double) ((int) unpacked[i]);
      }
    }
  }
  influx.close();
  cout<<"[done].\n";
  
  // PASS 2
  string outfile = outfix+".ldcor.pairs";
  ofstream fileOut(outfile.c_str());
  cout<<"# Pass 2: calculating LD correlations...";
  influx.open(bedfile.c_str(), std::ios::in | std::ios::binary);
  if(!influx){
    cerr << "[readGenotypes] Error reading file "<<bedfile<<endl;
    exit(1);
  }
  influx.seekg(0, ifstream::end);
  influx.seekg(3, ifstream::beg);
  double x;
  double sx, sx2, sy, sy2, sxy;
  double mx, mx2, my, my2, mxy;
  double num, den;
  double r, r_sq, r_sq_corr;
  double freqA1, freqa1;
  
  fileOut<<"CHR\tSNP1\tSNP2\tPOS1_KB\tPOS2_KB\tA1_1\tA1_2\tA2_1\tA2_2\tFREQ_A1_1\tFREQ_A1_2\tR_LD\tR_LD_SQ\tR_LD_SQ_CORR\n";
  for(j=0;j<M;j++){
    influx.read((char*)packed, sizeof(char) * numBytes);
    for(k=0;k<nSNPfound;k++){
      if(CHROM[j]==chrom[k] and fabs(POSKB[j]-poskb[k])<=window_kb){
        decode_plink(unpacked, packed, numBytes);
        nEff = 0;
        sxy  = 0.;
        sx2  = 0.;
        sy2  = 0.;
        sx   = 0.;
        sy   = 0.;
        for(i=0;i<N;i++){
          x = (double) ((int) unpacked[i]);
          if(x!=3. and xTarget[i+k*N]!=3){
            nEff++;
            sx  += xTarget[i+k*N];
            sy  += x;
            sx2 += xTarget[i+k*N]*xTarget[i+k*N];
            sy2 += x * x;
            sxy += x * xTarget[i+k*N];
          }
        }
        mx  = sx/nEff;
        my  = sy/nEff;
        mxy = sxy/nEff;
        mx2 = sx2/nEff;
        my2 = sy2/nEff;
        num = mxy-mx*my;
        den = sqrt((mx2-mx*mx)*(my2-my*my));
        
        freqa1 = 0.5*mx;
        freqA1 = 0.5*my;
        
        if(den > 0.f){
          r         = num / den;
          r_sq      = r * r ;
          r_sq_corr = r_sq - 1./(nEff-3);
          
          fileOut<<CHROM[j]<<"\t";
          fileOut<<snps[k]<<"\t";
          fileOut<<SNPS[j]<<"\t";
          fileOut<<poskb[k]<<"\t";
          fileOut<<POSKB[j]<<"\t";
          fileOut<<a1[k]<<"\t";
          fileOut<<A1[j]<<"\t";
          fileOut<<a2[k]<<"\t";
          fileOut<<A2[j]<<"\t";
          fileOut<<freqa1<<"\t";
          fileOut<<freqA1<<"\t";
          fileOut<<r<<"\t";
          fileOut<<r_sq<<"\t";
          fileOut<<r_sq_corr<<"\n";
        }
      }
    }
  }
  fileOut.close();
  cout<<"[done].\n";
  
  delete [] packed;
  delete [] unpacked;    
  delete [] CHROM;
  delete [] chrom;
  delete [] SNPS;
  delete [] snps;
  delete [] POSKB;
  delete [] poskb;
  delete [] A1;
  delete [] a1;
  delete [] A2;
  delete [] a2;
  delete [] snplist;
  delete [] snpsFound;
  delete [] xTarget;
  return EXIT_SUCCESS;
}


