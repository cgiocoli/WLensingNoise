#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <CCfits/CCfits>
#include "noisekappa.h"

using namespace CCfits;

void writeFits(std:: string filename,std:: valarray<float> f, int npix, int npixy,double zs){
  filename = "!"+filename;
  long naxis=2;
  long naxes[2]={npix,npixy};
  std::auto_ptr<FITS> fout(new FITS(filename,FLOAT_IMG,naxis,naxes));
  std::vector<long> naxex(2);
  naxex[0]=npix;
  naxex[1]=npixy;
  PHDU *phout=&fout->pHDU();
  phout->write( 1, npix*npixy, f );
  phout->addKey ("ZSOURCE",zs, "source redshift");
}

int main(){
  time_t start1,end1;
  double dif_time; 
  
  int n=2048;
  
  std:: valarray<float> km(n*n);
  
  double boxlrad=5*M_PI/180.;  // radiants

  double *noisekappa;
  noisekappa = new double[n*n]; // squared map!!!
  std:: cout << boxlrad << "  " << n << "  " << n << std:: endl;
  double zs;
  std:: vector<double> Pk, ll;
  std:: string fileink;
  std:: cin >> fileink;
  std:: cin >> zs;
  long seedi;
  std:: cin >> seedi;  
  std:: ifstream Kpl;  
  Kpl.open(fileink.c_str());
  if(Kpl.is_open()){
    /* build the noisy map if a file KappaPowerSpec.dat exists */
    double butr1,butr2;
    int nlines=0;
    while(Kpl >> butr1 >> butr2){
      Pk.push_back(butr2/4./pi/pi);
      ll.push_back(butr1);
      nlines++;
    }
    Kpl.close();
  }else{
    std:: cout << " file with Pkl does not exsists! " << std:: endl;
    exit(1);
  }

  noisekappa = getkappanoise(Pk,ll,boxlrad,n,seedi);  
  
  for(int j=0;j<n;j++){
    for(int i=0;i<n;i++){      
      km[i+n*j]+=float(noisekappa[i+n*j]);
    }
  }
  delete [] noisekappa;

  std:: string filout = "noise_map.fits";
  //writeFits(filout,km,n,n,zs);
  double *lll;
  double *Pl;
  int nb = 128;
  lll=new double[nb];
  Pl=new double[nb];  
  
  powerl(km,km,n,n,boxlrad,boxlrad,lll,Pl,nb);  
  std:: ofstream outfile;
  outfile.open(filout+"_mapPowerSpectrum.dat");
  for(int i=0;i<nb;i++){
    outfile << lll[i] << "  " << Pl[i] << std::endl;
  }
  outfile.close();    
}
