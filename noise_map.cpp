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

float getY(std:: vector<float> x, std:: vector<float> y,float xi){
  int nn = x.size();                                                               
  if(x[0]<x[nn-1]){                                                                                     
    if(xi>x[nn-1]) return y[nn-1];
    if(xi<x[0]) return y[0];
  }      
  else{
    if(xi<x[nn-1]) return y[nn-1];
    if(xi>x[0]) return y[0];        
  }                                    
  int i = locate (x,xi);           
  i = std::min (std::max (i,0), int (nn)-2);
  double f=(xi-x[i])/(x[i+1]-x[i]);
  if(i>1 && i<nn-2){
    double a0,a1,a2,a3,f2;                                                       
    f2 = f*f;                                                    
    a0 = y[i+2] - y[i+1] - y[i-1] + y[i];            
    a1 = y[i-1] - y[i] - a0;                                                                              
    a2 = y[i+1] - y[i-1];                                                              
    a3 = y[i];                                                                                     
    return a0*f*f2+a1*f2+a2*f+a3;                                                    
  }                                                                      
  else return f*y[i+1]+(1-f)*y[i];           
} 

void writeFits(std:: string filename,std:: valarray<float> f, int npix, int npixy,double zs,double fxrad, double fyrad,
	       double ra_c, double dec_c){
  // double fyrad=fxrad/npix*npixy;
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
  phout->addKey ("RADESYS" , "FK5" , "Coordinate system");
  phout->addKey ("CTYPE1" , "RA---TAN" , "Coordinates -- projection");
  phout->addKey ("CTYPE2" , "DEC--TAN" , "Coordinates -- projection");
  phout->addKey ("EQUINOX" , 2000 , "Epoch of the equinox");
  phout->addKey ("CUNIT1" , "deg" , " ");
  phout->addKey ("CUNIT2" , "deg" , " ");
  // if even
  double xo = double(npix)/2.+1;
  double yo = double(npixy)/2.+1;
  phout->addKey ("CRPIX1",xo,"X reference pixel");
  phout->addKey ("CRPIX2",yo,"Y reference pixel");
  double xov= ra_c;
  double yov= dec_c;
  phout->addKey ("CRVAL1",xov,"Reference longitude");
  phout->addKey ("CRVAL2",yov,"Reference latitude");
  double dltx = -fxrad/npix*180./M_PI;
  double dlty = fyrad/npixy*180./M_PI;
  phout->addKey ("CDELT1",dltx,"X scale");
  phout->addKey ("CDELT2",dlty,"Y scale");  
}

void readFits (std::string fn, valarray<float> &map, int &nx, int &ny, double &ra_c, double &dec_c){
  std::auto_ptr<FITS> ff(new FITS (fn, Read));
  PHDU *h0=&ff->pHDU();
  // number of pixels x and y (usually MOKA are squares)
  nx=h0->axis(0);
  ny=h0->axis(1);
  // read the fits file
  h0->read(map);
  try {
    h0->readKey ("CRVAL1",ra_c);
  }
  catch(CCfits::HDU::NoSuchKeyword){
    ra_c=0;
  }  
  try {
    h0->readKey ("CRVAL2",dec_c);
  }
  catch(CCfits::HDU::NoSuchKeyword){
    dec_c=0;
  }  
}

string getFileName(const string& s) {
  char sep = '/';
  
#ifdef _WIN32
  sep = '\\';
#endif
  size_t i = s.rfind(sep, s.length());
  if (i != string::npos) {
    return(s.substr(i+1, s.length() - i));
  }
  
  return("");
}

int main(){
  time_t start1,end1;
  double dif_time;
  
  double range;
  std:: cin >> range;
  double boxlrad=range*M_PI/180.;  // radiants

  double boxlradx = boxlrad;
  double boxlrady = boxlradx;  // squared map!!!
  double zs;
  std:: vector<double> Pk, ll;
  std:: string fileink;
  std:: cin >> fileink;
  std:: cin >> zs;
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
  std:: string filin;
  std:: cin >> filin;
  std:: string filplmap;
  std:: cin >> filplmap;  
  long seedi;
  std:: cin >> seedi;
  std:: valarray<float> kin;
  double filter;
  std:: cin >> filter;
  std:: string filin0 = getFileName(filin);
  if(filin0=="") filin0=filin;
  double ra_c,dec_c;
  int n;  
  readFits(filin,kin,n,n,ra_c,dec_c);
  
  std:: valarray<float> km(n*n);  

  double *noisykappa;
  noisykappa = new double[n*n]; // squared map!!!
  std:: cout << boxlrad << "  " << n << "  " << n << std:: endl;
  
  // reading the power spectrum
  std:: vector<double> Pkm, llm;  
  Kpl.open(filplmap.c_str());
  if(Kpl.is_open()){
    /* build the noisy map if a file KappaPowerSpec.dat exists */
    double butr1,butr2;
    int nlines=0;
    while(Kpl >> butr1 >> butr2){
      Pkm.push_back(butr2);
      llm.push_back(butr1);
      nlines++;
    }
    Kpl.close();
  }else{
    std:: cout << " file with Pkl for the map does not exsists! " << std:: endl;
    exit(1);
  }
  if(filter<-0.5){
    noisykappa = getkappanoise(Pk,ll,boxlrad,n,seedi);    
  }else{
    noisykappa = getkappanoise3(Pk,ll,boxlrad,n,kin,Pkm,llm,filter,seedi);
  }
  
  for(int j=0;j<n;j++){
    for(int i=0;i<n;i++){      
      km[i+n*j]+=float(noisykappa[i+n*j]);
    }
  }
  delete [] noisykappa;
  std:: string filout = "noise_map.fits4"+filin0;  
  writeFits(filout,km,n,n,zs,boxlradx,boxlrady,ra_c,dec_c);
  double *lll;
  double *Pl;
  int nb = 128;
  lll=new double[nb];
  Pl=new double[nb];  
  
  powerl(km,km,n,n,boxlradx,boxlrady,lll,Pl,nb);  
  std:: ofstream outfile;
  outfile.open(filout+"_mapPowerSpectrum.dat");
  for(int i=0;i<nb;i++){
    outfile << lll[i] << "  " << Pl[i] << std::endl;
  }
  outfile.close();  

  for(int j=0;j<n;j++){
    for(int i=0;i<n;i++){      
      km[i+n*j]=km[i+n*j]+kin[i+n*j];
    }
  }

  float averagek = km.sum()/float(n)/float(n);
  for(int j=0;j<n;j++){
    for(int i=0;i<n;i++){      
      km[i+n*j]=km[i+n*j]-averagek;
    }
  } 
  if(filter<-0.5){
    filout = filin0+"_noised_uncorrelated.fits";
  }else{
    filout = filin0+"_noised.fits";
  }
  writeFits(filout,km,n,n,zs,boxlradx,boxlrady,ra_c,dec_c);
  powerl(km,km,n,n,boxlradx,boxlrady,lll,Pl,nb);  
  // std:: ofstream outfile;
  if(filter<-0.5){
    outfile.open(filin0+"_noised_uncorrelated_mapPowerSpectrum.dat");
  }else{
    outfile.open(filin0+"_noised_mapPowerSpectrum.dat");
  }
  for(int i=0;i<nb;i++){
    outfile << lll[i] << "  " << Pl[i] << std::endl;
  } 
  outfile.close();
  return 0;  
}



