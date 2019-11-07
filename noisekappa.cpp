#include "noisekappa.h"
#include <time.h>
#include <random>

int MULTIPLIER=3;
int DIVISOR=7;




double* getkappanoise3(std:: vector<double> Pklin, std:: vector<double> llin,
		       double boxlrad,int npix, std:: valarray<float> map,
		       std:: vector<double> Pk, std:: vector<double> l,
		       double filter,long seedi){
  
  double sigmag = filter*M_PI/180.;  
  // generate a map for the noise
  double *input=new double[npix*npix];
  fftw_complex *output=new fftw_complex[npix*(npix/2+1)];
  for (int i=0; i<npix*npix; i++) input[i] = double(map[i]);
  fftw_plan pforward;  
  pforward=fftw_plan_dft_r2c_2d(npix,npix,input,output,FFTW_ESTIMATE);
  fftw_execute( pforward );
  // this part is here to initialize the output map
  
  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  double mean = 0.0;
  double stddev  = 1.0;
  std::normal_distribution<double> normal(mean, stddev);

  long seed = time(NULL) + MULTIPLIER * clock() % DIVISOR;
  if(seedi>0){
    seed = seedi;
    gsl_rng_set(r,seed);
  }else{
    gsl_rng_set(r,34872); // <--- if you want to set the seed by hand
  }
  int i,j,ii;
  double k_dir[2];
  double lmag, p;
  fftw_plan pbackward;

  double value, phase;
  fftw_complex *kappak;
  kappak = new fftw_complex[npix*(npix/2+1)];
  double K = 2. * M_PI  / boxlrad;    
  
  for (i=0; i<npix; i++) for (j=0; j<npix/2+1; j++){
      if(i < npix / 2) k_dir[0] = i*K;
      else k_dir[0] = -(npix-i)*K;
      if(j < npix / 2) k_dir[1] = j*K;
      else k_dir[1] = -(npix-j)*K;
      ii = i*(npix/2+1) + j;
      
      lmag = sqrt(k_dir[0]*k_dir[0] + k_dir[1]*k_dir[1]);	
      if((j == 0 && i == 0) || lmag / K > npix / 2.){
	kappak[ii][0] = 0.0;
	kappak[ii][1] = 0.0;
      }
      else{
	p = getY(llin,Pklin,lmag);       
 	// double gg = (normal(generator));
	value = fabs(gsl_ran_ugaussian(r)) * sqrt(p) * K;
	//value = gg * sqrt(p) * K;		
	phase = 2.0 * M_PI * gsl_rng_uniform(r);
	double fase0,fase1;
	// this is different from 0 if you want to create a map in phase with the input
	double val0 = sqrt(output[ii][0]*output[ii][0] + output[ii][1]*output[ii][1]);
	if(val0>1e-170){
	  // in phase
	  fase0 = output[ii][0]/val0;
	  fase1 = output[ii][1]/val0;
	  phase = atan2(fase1,fase0);
	  fase0 = cos(phase);
	  fase1 = sin(phase);
	}else{
	  // random
	  fase0 = cos(phase);
	  fase1 = sin(phase);		  
	}
	// set the value
	kappak[ii][0] = value * fase0;
	kappak[ii][1] = value * fase1;		

	if(fase0!=fase0){
	  std:: cout << fase0 << "  " << fase1 << std:: endl;	  
	  exit(1);
	}
	if(fase1!=fase1){
	  std:: cout << fase0 << "  " << fase1 << std:: endl;	  
	  exit(1);
	}
	if(kappak[ii][0]!=kappak[ii][0]){
	  std:: cout << kappak[ii][0] << std:: endl;
	  exit(1);
	}
	if(kappak[ii][1]!=kappak[ii][1]){
	  std:: cout << kappak[ii][1] << std:: endl;
	  exit(1);
	}
      }
    }
  
  double *kappa=new double[npix*npix];      
  pbackward = fftw_plan_dft_c2r_2d(npix, npix, kappak, kappa, FFTW_ESTIMATE);
  fftw_execute(pbackward);
  
  // now constract the sum of the maps
  std:: valarray<float> mapt(npix*npix);
  std:: valarray<float> mapkappa(npix*npix);  
  for(int j=0;j<npix;j++) for(int i=0;i<npix;i++) {
      mapt[i+npix*j] = map[i+npix*j] + kappa[i+npix*j];
      mapkappa[i+npix*j] = kappa[i+npix*j];
    }  
  // compute the power spectrum of the sum
  double *lt;
  double *Plt;
  int nb = l.size();
  lt=new double[nb];
  Plt=new double[nb];   
  powerl(mapt,mapt,npix,npix,boxlrad,boxlrad,lt,Plt,nb);
  // compute the power spectrum of the gaussian linear map
  double *l0;
  double *Pl0;
  l0=new double[nb];
  Pl0=new double[nb];  
  powerl(mapkappa,mapkappa,npix,npix,boxlrad,boxlrad,l0,Pl0,nb);
  
  std:: vector<double> A,lA;
  for(int i=0;i<nb;i++){
    if(Plt[i]>0){
      A.push_back(Pl0[i]/Plt[i]);
      lA.push_back(l[i]);    
      std:: cout << " amplitude " <<
	Plt[i] << "  " << Pl0[i] << "  " << Pk[i] << "  "
		 << l[i] << "  " << Pl0[i]/Plt[i] << std:: endl;
    }
  }

  std:: cout << " done .... map " << std:: endl;

  gsl_rng *rr = gsl_rng_alloc(gsl_rng_mt19937);
  // now build the map considering the factor
  if(seedi>0){
    seed = seedi;
    gsl_rng_set(rr,seed);     
  }else{
    gsl_rng_set(rr,34872); // <--- if you want to set the seed by hand
  }  
  std:: cout << " seed reset " << std:: endl;
  
  for (i=0; i<npix; i++) for (j=0; j<npix/2+1; j++){
      if(i < npix / 2) k_dir[0] = i*K;
      else k_dir[0] = -(npix-i)*K;
      if(j < npix / 2) k_dir[1] = j*K;
      else k_dir[1] = -(npix-j)*K;
      ii = i*(npix/2+1) + j;
      
      lmag = sqrt(k_dir[0]*k_dir[0] + k_dir[1]*k_dir[1]);	
      if((j == 0 && i == 0) || lmag / K > npix / 2.){
	kappak[ii][0] = 0.0;
	kappak[ii][1] = 0.0;
      }
      else{
	p = getY(llin,Pklin,lmag);
	double pfactor = getY(lA,A,lmag);
	//double gg = normal(generator2);
	value = fabs(gsl_ran_ugaussian(rr)) * sqrt(p*pfactor) * K;
	//value = gg * sqrt(p*pfactor) * K;	
	phase = 2.0 * M_PI * gsl_rng_uniform(rr);
	double fase0,fase1;
	// this is different from 0 if you want to create a map in phase with the input	
	double val0 = sqrt(output[ii][0]*output[ii][0] + output[ii][1]*output[ii][1]);
	if(val0>1e-170){
	  // in phase
	  fase0 = output[ii][0]/val0;
	  fase1 = output[ii][1]/val0;
	  phase = atan2(fase1,fase0);
	  fase0 = cos(phase);
	  fase1 = sin(phase);
	}else{
	  // random
	  fase0 = cos(phase);
	  fase1 = sin(phase);		  
	}
	// apply smoothing if larger than 0
	kappak[ii][0] = value * fase0 * exp(-2*M_PI*M_PI*lmag*lmag*sigmag*sigmag);
	kappak[ii][1] = value * fase1 * exp(-2*M_PI*M_PI*lmag*lmag*sigmag*sigmag);
	
	if(fase0!=fase0){
	  std:: cout << fase0 << "  " << fase1 << std:: endl;	  
	  exit(1);
	}
	if(fase1!=fase1){
	  std:: cout << fase0 << "  " << fase1 << std:: endl;	  
	  exit(1);
	}
	if(kappak[ii][0]!=kappak[ii][0]){
	  std:: cout << "1 " << kappak[ii][0] << std:: endl;
	  exit(1);
	}
	if(kappak[ii][1]!=kappak[ii][1]){
	  std:: cout << "2 " << kappak[ii][1] << std:: endl;
	  exit(1);
	}
      }
    }
  
  pbackward = fftw_plan_dft_c2r_2d(npix, npix, kappak, kappa, FFTW_ESTIMATE);
  fftw_execute(pbackward);
  fftw_destroy_plan(pbackward);
  gsl_rng_free (rr);        
  gsl_rng_free (r);        
  fftw_free(kappak);
  
  delete[] input;
  delete[] output;
  fftw_destroy_plan(pforward);
  return kappa;
}

double* getkappanoise(std:: vector<double> Pklin, std:: vector<double> llin,
		      double boxlrad,int npix,long seedi){


  // generate a random map for the noise
  fftw_complex *output=new fftw_complex[npix*(npix/2+1)];
  // this part initialize the output map all zero
  
  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);  
  double mean = 0.0;
  double stddev  = 1.0;
  std::normal_distribution<double> normal(mean, stddev);
  
  long seed = time(NULL) + MULTIPLIER * clock() % DIVISOR;
  if(seedi>0){
    seed = seedi;
    gsl_rng_set(r,seed);     
  }else{
    gsl_rng_set(r,34872); // <--- if you want to set the seed by hand
  }
  int i,j,ii;
  double k_dir[2];
  double lmag, p;
  fftw_plan pbackward;

  double value, phase;
  fftw_complex *kappak;
  kappak = new fftw_complex[npix*(npix/2+1)];
  double K = 2. * M_PI  / boxlrad;    
  
  for (i=0; i<npix; i++) for (j=0; j<npix/2+1; j++){
      if(i < npix / 2) k_dir[0] = i*K;
      else k_dir[0] = -(npix-i)*K;
      if(j < npix / 2) k_dir[1] = j*K;
      else k_dir[1] = -(npix-j)*K;
      ii = i*(npix/2+1) + j;
      
      lmag = sqrt(k_dir[0]*k_dir[0] + k_dir[1]*k_dir[1]);	
      if((j == 0 && i == 0) || lmag / K > npix / 2.){
	kappak[ii][0] = 0.0;
	kappak[ii][1] = 0.0;
      }
      else{
	p = getY(llin,Pklin,lmag);
	// double gg = normal(generator);
	value = fabs(gsl_ran_ugaussian(r)) * sqrt(p) * K;
	// value = fabs(gg) * sqrt(p) * K;	
	phase = 2.0 * M_PI * gsl_rng_uniform(r);
	double fase0,fase1;
	// in this cae val0 is always 0	
	double val0 = sqrt(output[ii][0]*output[ii][0] + output[ii][1]*output[ii][1]);
	if(val0>1e-170){
	  //  in phase
	  fase0 = output[ii][0]/val0;
	  fase1 = output[ii][1]/val0;
	  phase = atan2(fase1,fase0);
	  fase0 = cos(phase);
	  fase1 = sin(phase);
	}else{
	  // random
	  fase0 = cos(phase);
	  fase1 = sin(phase);		  
	}	
	kappak[ii][0] = value * fase0;
	kappak[ii][1] = value * fase1;	
	
	if(fase0!=fase0){
	  std:: cout << fase0 << "  " << fase1 << std:: endl;	  
	  exit(1);
	}
	if(fase1!=fase1){
	  std:: cout << fase0 << "  " << fase1 << std:: endl;	  
	  exit(1);
	}
	if(kappak[ii][0]!=kappak[ii][0]){
	  std:: cout << kappak[ii][0] << std:: endl;
	  exit(1);
	}
	if(kappak[ii][1]!=kappak[ii][1]){
	  std:: cout << kappak[ii][1] << std:: endl;
	  exit(1);
	}
      }
    }
  
  double *kappa=new double[npix*npix];      
  pbackward = fftw_plan_dft_c2r_2d(npix, npix, kappak, kappa, FFTW_ESTIMATE);
  fftw_execute(pbackward);

  fftw_destroy_plan(pbackward);
  gsl_rng_free (r);        
  fftw_free(kappak);
  
  delete[] output;
  return kappa;
}
