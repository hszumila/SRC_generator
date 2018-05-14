 /* generates lund dat file for input to gemc              *
  * assumes scattering off 12C                             *
  * beamP is beam momentum, zVtx is vtx posn in cm,        *
  * n is number of events                                  *
  * run this in root as .x generator.C(beamP, zVtx, n)     */

#include <iostream>
#include "time.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"

const double PI = 3.14159265358979323846264338;
const double d2r = PI/180.0;
const double r2d = 180.0/PI;
const double MP = 0.938;
const double ME = 0.000511;
const double MA = 12.0*0.9315; //mass of carbon nucleus
const double MAM2 = 10.0*0.9315; //mass of carbon nucleus after losing a SRC pair 

Double_t fFormFactor(double eprime, double beamE, double theta){

  double a2 = 0.71;//GeV^2
  double mu_p = 2.79;
  
  theta = theta*d2r;
  double q2 = 4.0*eprime*beamE*pow(sin(theta/2.0),2.0);
  double xb = q2/(2.0*MP*(beamE-eprime));
  double Ge = pow(1.0+q2/a2,-2.0);
  double Gm = mu_p*Ge;
  double tau = q2/pow(2.0*MP,2.0);

  Double_t ff = (eprime/beamE)*((pow(Ge,2.0)+tau*pow(Gm,2.0))/(1.0+tau) + 2*tau*pow(Gm*tan(theta/2.0),2.0));
  return ff;
}
 
Double_t fMott(double theta, double beamE){

  Double_t ff = pow((1./137.0)*cos((theta*d2r)/2.0),2.0)/(4*pow(beamE,2.0)*pow(sin((theta*d2r)/2.0),4.0));
  return ff;  

}

//input magnitude of pmiss momentum and emiss in GeV
Double_t calcPrecoil(double pmiss, double emiss){
  double beta = MA - emiss; 
  double betaTilde = -0.5*(pow(pmiss,2.0)+pow(MAM2,2.0)-pow(MP,2.0)-pow(beta,2.0))/beta;
  double sol1 = (2*betaTilde*pmiss/beta + sqrt(4*pow(betaTilde/beta,2.0) - 4.0*(1.0-pow(pmiss/beta,2.0))*(pow(MP,2.0)-pow(betaTilde,2.0))))/(2.0*(1.0-pow(pmiss/beta,2.0)));
  double sol2 = (2*betaTilde*pmiss/beta - sqrt(4*pow(betaTilde/beta,2.0) - 4.0*(1.0-pow(pmiss/beta,2.0))*(pow(MP,2.0)-pow(betaTilde,2.0))))/(2.0*(1.0-pow(pmiss/beta,2.0)));

  return TMath::Max(sol1,sol2);

}

void generatorSRC(double beamP, double zVtx, int n)//input number of events
{
  srand (time(NULL)); //intialize random seed
  //FILE *outLund = fopen(Form("src_%.1fgev_%.1fcm_src.dat",beamE,zVtx),"a+");
  FILE *outLund = fopen("input.dat","a+");
  TFile *fout = new TFile(Form("src_%.1fgev_%.1fcm_src.root",beamP,zVtx),"RECREATE");
  TTree tt("tt","variable tree");

  double xB, Q2, ffWt, mottWt, theta_eprime, phi_eprime, p_eprime;
  double theta_pmiss, phi_pmiss, p_pmiss, km4Wt, e_pmiss;
  double theta_pcm, phi_pcm, p_pcm, e_pcm;
  double theta_precoil, phi_precoil, p_precoil, e_precoil;
  double theta_q, phi_q, p_q, omega;
  double theta_proton, phi_proton, p_proton, e_proton;
  double theta_p_q, theta_recoil_q, theta_miss_q, theta_miss_recoil, ratio_p_q;
  double mmiss;
  double epvec[3], pmvec[3], cmvec[3];

  tt.Branch("xB",&xB,"xB/D");
  tt.Branch("Q2",&Q2,"Q2/D");
  tt.Branch("ffWt",&ffWt,"ffWt/D");
  tt.Branch("mottWt",&mottWt,"mottWt/D");
  tt.Branch("theta_q",&theta_q,"theta_q/D");
  tt.Branch("phi_q",&phi_q,"phi_q/D");
  tt.Branch("p_q",&p_q,"p_q/D");
  tt.Branch("omega",&omega,"omega/D");
  tt.Branch("theta_eprime",&theta_eprime,"theta_eprime/D");
  tt.Branch("phi_eprime",&phi_eprime,"phi_eprime/D");
  tt.Branch("p_eprime",&p_eprime,"p_eprime/D");
  tt.Branch("theta_pmiss",&theta_pmiss,"theta_pmiss/D");
  tt.Branch("phi_pmiss",&phi_pmiss,"phi_pmiss/D");
  tt.Branch("p_pmiss",&p_pmiss,"p_pmiss/D");
  tt.Branch("e_pmiss",&e_pmiss,"e_pmiss/D");
  tt.Branch("km4Wt",&km4Wt,"km4Wt/D");
  tt.Branch("theta_pcm",&theta_pcm,"theta_pcm/D");
  tt.Branch("phi_pcm",&phi_pcm,"phi_pcm/D");
  tt.Branch("p_pcm",&p_pcm,"p_pcm/D");
  tt.Branch("e_pcm",&e_pcm,"e_pcm/D");
  tt.Branch("theta_precoil",&theta_precoil,"theta_precoil/D");
  tt.Branch("phi_precoil",&phi_precoil,"phi_precoil/D");
  tt.Branch("p_precoil",&p_precoil,"p_precoil/D");
  tt.Branch("e_precoil",&e_precoil,"e_precoil/D");
  tt.Branch("theta_proton",&theta_proton,"theta_proton/D");
  tt.Branch("phi_proton",&phi_proton,"phi_proton/D");
  tt.Branch("p_proton",&p_proton,"p_proton/D");
  tt.Branch("e_proton",&e_proton,"e_proton/D");
  tt.Branch("theta_p_q",&theta_p_q,"theta_p_q/D");
  tt.Branch("theta_recoil_q",&theta_recoil_q,"theta_recoil_q/D");
  tt.Branch("theta_miss_q",&theta_miss_q,"theta_miss_q/D");
  tt.Branch("theta_miss_recoil",&theta_miss_recoil,"theta_miss_recoil/D");
  tt.Branch("ratio_p_q",&ratio_p_q,"ratio_p_q/D");
  tt.Branch("mmiss",&mmiss,"mmiss/D");
  
  //////////////////
  //write the beamE
  //////////////////
  TLorentzVector eBeam(0.0, 0.0, beamP, sqrt(pow(beamP,2.0)+pow(ME,2.0)));

  for (int i=0; i<n; i++){

    //////////////////////////////////////////////////////////////////
    //generate E', calc Mott weight, calc FF weight, calc Q2, calc xB
    /////////////////////////////////////////////////////////////////
    p_eprime = gRandom->Uniform(0.1,beamP);
    if (i==0){cout<<"First event: "<<p_eprime<<" GeV."<<endl;}
    gRandom->Sphere(epvec[0],epvec[1],epvec[2],p_eprime);
    TLorentzVector eprime(epvec[0], epvec[1], epvec[2], sqrt(pow(p_eprime,2.0)+pow(ME,2.0)));
    theta_eprime = acos(eprime.Pz()/eprime.P())*r2d;
    phi_eprime = atan(eprime.Py()/eprime.Px())*r2d;
  
    TLorentzVector qvec = eBeam - eprime;
    omega = qvec.E();
    theta_q = acos(qvec.Pz()/qvec.P())*r2d;
    phi_q = atan(qvec.Py()/qvec.Px())*r2d;
    p_q = qvec.P();
    
    mottWt = fMott(theta_eprime, beamP);
    ffWt = fFormFactor(p_eprime, beamP, theta_eprime);
    Q2 = 4*beamP*p_eprime*pow(sin(theta_eprime*d2r/2.0),2.0);
    xB = Q2/(2.0*MP*omega);
    
    ///////////////////////////////////
    //generate Pmiss, calc 1/k4 weight
    ///////////////////////////////////
    p_pmiss = gRandom->Uniform(0.25,1.0);
    gRandom->Sphere(pmvec[0],pmvec[1],pmvec[2],p_pmiss);
    TLorentzVector pmiss(pmvec[0], pmvec[1], pmvec[2], 0.0);
    theta_pmiss = acos(pmiss.Pz()/pmiss.P())*r2d;
    phi_pmiss = atan(pmiss.Py()/pmiss.Px())*r2d;
    //km4Wt = (p_pmiss<0.25)?153.6*pow(p_pmiss,2.0):0.0666666667/pow(p_pmiss,4.0);
    km4Wt = 0.0666666667/pow(p_pmiss,4.0);

    ////////////////////////////////
    //calculate the proton and Emiss
    ////////////////////////////////
    TLorentzVector proton = qvec + pmiss;
    e_proton = sqrt(pow(proton.P(),2.0) + pow(MP,2.0));
    theta_proton = acos(proton.Pz()/proton.P())*r2d;
    phi_proton = atan(proton.Py()/proton.Px())*r2d;
    p_proton = proton.P();
    proton.SetE(e_proton);
    e_pmiss = proton.E() - qvec.E();
    pmiss.SetE(e_pmiss);
    
    ////////////////////////////////////
    //calculate the magnitude of Precoil
    ////////////////////////////////////
    p_precoil = calcPrecoil(pmiss.P(),pmiss.E());
    p_pcm = sqrt(pow(p_pmiss,2.0)+pow(p_precoil,2.0)-2.0*p_pmiss*p_precoil);
   
    ///////////////
    //generate Pcm
    ///////////////
    gRandom->Sphere(cmvec[0],cmvec[1],cmvec[2],p_pcm);
    e_pcm = MA-sqrt(pow(p_pcm,2.0)+pow(MAM2,2.0));
    TLorentzVector pcm(cmvec[0], cmvec[1], cmvec[2], e_pcm);
    theta_pcm = acos(pcm.Pz()/pcm.P())*r2d;
    phi_pcm = atan(pcm.Py()/pcm.Px())*r2d;
    
    ////////////////////////////////
    //calculate Precoil and Erecoil
    ////////////////////////////////
    TLorentzVector precoil = pcm - pmiss;
    e_precoil = precoil.E();
    theta_precoil = acos(precoil.Pz()/precoil.P())*r2d;
    phi_precoil = atan(precoil.Py()/precoil.Px())*r2d;

    ////////////////////////////////
    // calc other useful quantities
    ////////////////////////////////
    theta_p_q = acos((qvec.Px()*proton.Px()+qvec.Py()*proton.Py()+qvec.Pz()*proton.Pz())/(qvec.P()*proton.P()))*r2d;
    theta_recoil_q = acos((qvec.Px()*precoil.Px()+qvec.Py()*precoil.Py()+qvec.Pz()*precoil.Pz())/(qvec.P()*precoil.P()))*r2d;
    theta_miss_q = acos((qvec.Px()*pmiss.Px()+qvec.Py()*pmiss.Py()+qvec.Pz()*pmiss.Pz())/(qvec.P()*pmiss.P()))*r2d;
    theta_miss_recoil = acos((pmiss.Px()*precoil.Px()+pmiss.Py()*precoil.Py()+pmiss.Pz()*precoil.Pz())/(pmiss.P()*precoil.P()))*r2d;
    ratio_p_q = (proton.P()/qvec.P());

    //TLorentzVector deuteron(0,0,0,1.8756);
    mmiss = sqrt(pow(pmiss.E(),2.0) - pow(pmiss.P(),2.0));
    
    if (ratio_p_q>0.6 && ratio_p_q<1.0 && xB>0.1 &&xB<2.0 && p_pcm<0.25 && theta_eprime<45.0){
      if (i%1000==0){cout<<i<<" events"<<endl;}
      tt.Fill();
      
      fprintf(outLund, "3\t 0.\t 0.\t 0.\t 0.5\t 11.\t %f\t 0.\t 0.\t 1.\n", beamP);
      fprintf(outLund, "1\t -1.\t 1.\t 11.\t 0.\t 0.\t %f\t %f\t %f\t %f\t 0.000511\t 0.0\t 0.0\t %f\n",eprime.Px(),eprime.Py(),eprime.Pz(),eprime.E(),zVtx);
      fprintf(outLund, "2\t 1.\t 1.\t 2212.\t 0.\t 0.\t %f\t %f\t %f\t %f\t 0.9383\t 0.0\t 0.0\t %f\n",proton.Px(),proton.Py(),proton.Pz(),proton.E(),zVtx);
      fprintf(outLund, "3\t 1.\t 1.\t 2212.\t 0.\t 0.\t %f\t %f\t %f\t %f\t 0.9383\t 0.0\t 0.0\t %f\n",precoil.Px(),precoil.Py(),precoil.Pz(),precoil.E(),zVtx);
      
    }
  else {i--;}

  }

  tt.Write();   
  fclose(outLund);
  
  fout->Write();
  fout->Close();

}
