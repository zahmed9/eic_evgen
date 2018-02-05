#include <iostream>
#include <iostream>     // std::cout
#include <fstream>      // std::ifstream
#include <stdio.h>
#include <cstdlib>
#include <stdlib.h>
#include <iomanip>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <string>       // std::string
#include <sstream>      // std::stringstream, std::stringbuf
#include "TStopwatch.h"
#include "TDatime.h"

void pim::setrootfile( string rootFile ){ 


///****************************************
/// Bill: re-delcreation of f1 is fixed,
///       all object function calls are switched pointer function calls


  f = new TFile(rootFile.c_str(),"recreate"); 

  t1 = new TTree();
  t1->SetName("t1");

  // 6 Particles, electron, proton, scat. electron, photon, pion and neutron in col frame
  //-----------------------------------------------------------------------------------------------------

  /* t1->Branch("Proton_Theta_Col",                          &fProton_Theta_Col,                          "fProton_Theta_Col/D"); */
  /* t1->Branch("Proton_Phi_Col",                            &fProton_Phi_Col,                            "fProton_Phi_Col/D"); */
  /* t1->Branch("Proton_Energy_Col_GeV",                     &fProton_Energy_Col_GeV,                     "fProton_Energy_Col_GeV/D"); */
  /* t1->Branch("Proton_Mom_Col_GeV",                        &fProton_Mom_Col_GeV,                        "fProton_Mom_Col_GeV/D"); */
  /* t1->Branch("Proton_MomZ_Col_GeV",                       &fProton_MomZ_Col_GeV,                       "fProton_MomZ_Col_GeV/D"); */
  /* t1->Branch("Proton_MomX_Col_GeV",                       &fProton_MomX_Col_GeV,                       "fProton_MomX_Col_GeV/D"); */
  /* t1->Branch("Proton_MomY_Col_GeV",                       &fProton_MomY_Col_GeV,                       "fProton_MomY_Col_GeV/D"); */

  /* t1->Branch("Electron_Theta_Col",                        &fElectron_Theta_Col,                        "fElectron_Theta_Col/D"); */
  /* t1->Branch("Electron_Phi_Col",                          &fElectron_Phi_Col,                          "fElectron_Phi_Col/D"); */
  /* t1->Branch("Electron_Energy_Col_GeV",                   &fElectron_Energy_Col_GeV,                   "fElectron_Energy_Col_GeV/D"); */
  /* t1->Branch("Electron_Mom_Col_GeV",                      &fElectron_Mom_Col_GeV,                      "fElectron_Mom_Col_GeV/D"); */
  /* t1->Branch("Electron_MomX_Col_GeV",                     &fElectron_MomX_Col_GeV,                     "fElectron_MomX_Col_GeV/D"); */
  /* t1->Branch("Electron_MomY_Col_GeV",                     &fElectron_MomY_Col_GeV,                     "fElectron_MomY_Col_GeV/D"); */
  /* t1->Branch("Electron_MomZ_Col_GeV",                     &fElectron_MomZ_Col_GeV,                     "fElectron_MomZ_Col_GeV/D"); */

  t1->Branch("ScatElec_Theta_Col",                        &fScatElec_Theta_Col,                        "fScatElec_Theta_Col/D");
  t1->Branch("ScatElec_Phi_Col",                          &fScatElec_Phi_Col,                          "fScatElec_Phi_Col/D");
  t1->Branch("ScatElec_Energy_Col_GeV",                   &fScatElec_Energy_Col_GeV,                   "fScatElec_Energy_Col_GeV/D");
  t1->Branch("ScatElec_Mom_Col_GeV",                      &fScatElec_Mom_Col_GeV,                      "fScatElec_Mom_Col_GeV/D");
  t1->Branch("ScatElec_MomX_Col_GeV",                     &fScatElec_MomX_Col_GeV,                     "fScatElec_MomX_Col_GeV/D");
  t1->Branch("ScatElec_MomY_Col_GeV",                     &fScatElec_MomY_Col_GeV,                     "fScatElec_MomY_Col_GeV/D");
  t1->Branch("ScatElec_MomZ_Col_GeV",                     &fScatElec_MomZ_Col_GeV,                     "fScatElec_MomZ_Col_GeV/D");

  /* t1->Branch("Photon_Theta_Col",                          &fPhoton_Theta_Col,                          "fPhoton_Theta_Col/D"); */
  /* t1->Branch("Photon_Phi_Col",                            &fPhoton_Phi_Col,                            "fPhoton_Phi_Col/D"); */
  /* t1->Branch("Photon_Energy_Col_GeV",                     &fPhoton_Energy_Col_GeV,                     "fPhoton_Energy_Col_GeV/D"); */
  /* t1->Branch("Photon_Mom_Col_GeV",                        &fPhoton_Mom_Col_GeV,                        "fPhoton_Mom_Col_GeV/D"); */
  /* t1->Branch("Photon_MomX_Col_GeV",                       &fPhoton_MomX_Col_GeV,                       "fPhoton_MomX_Col_GeV/D"); */
  /* t1->Branch("Photon_MomY_Col_GeV",                       &fPhoton_MomY_Col_GeV,                       "fPhoton_MomY_Col_GeV/D"); */
  /* t1->Branch("Photon_MomZ_Col_GeV",                       &fPhoton_MomZ_Col_GeV,                       "fPhoton_MomZ_Col_GeV/D"); */

  t1->Branch("Pion_Theta_Col",                            &fPion_Theta_Col,                            "fPion_Theta_Col/D");
  t1->Branch("Pion_Phi_Col",                              &fPion_Phi_Col,                              "fPion_Phi_Col/D");
  t1->Branch("Pion_Energy_Col_GeV",                       &fPion_Energy_Col_GeV,                       "fPion_Energy_Col_GeV/D");
  t1->Branch("Pion_Mom_Col_GeV",                          &fPion_Mom_Col_GeV,                          "fPion_Mom_Col_GeV/D");
  t1->Branch("Pion_MomX_Col_GeV",                         &fPion_MomX_Col_GeV,                         "fPion_MomX_Col_GeV/D");
  t1->Branch("Pion_MomY_Col_GeV",                         &fPion_MomY_Col_GeV,                         "fPion_MomY_Col_GeV/D");
  t1->Branch("Pion_MomZ_Col_GeV",                         &fPion_MomZ_Col_GeV,                         "fPion_MomZ_Col_GeV/D");

  t1->Branch("Neutron_Theta_Col",                         &fNeutron_Theta_Col,                         "fNeutron_Theta_Col/D");
  t1->Branch("Neutron_Phi_Col",                           &fNeutron_Phi_Col,                           "fNeutron_Phi_Col/D");
  t1->Branch("Neutron_Energy_Col_GeV",                    &fNeutron_Energy_Col_GeV,                    "fNeutron_Energy_Col_GeV/D");
  t1->Branch("Neutron_Mom_Col_GeV",                       &fNeutron_Mom_Col_GeV,                       "fNeutron_Mom_Col_GeV/D");
  t1->Branch("Neutron_MomX_Col_GeV",                      &fNeutron_MomX_Col_GeV,                      "fNeutron_MomX_Col_GeV/D");
  t1->Branch("Neutron_MomY_Col_GeV",                      &fNeutron_MomY_Col_GeV,                      "fNeutron_MomY_Col_GeV/D");
  t1->Branch("Neutron_MomZ_Col_GeV",                      &fNeutron_MomZ_Col_GeV,                      "fNeutron_MomZ_Col_GeV/D");

  // 6 Particles, electron, proton, scat. electron, photon, pion and neutron in proton's rest frame
  //-----------------------------------------------------------------------------------------------

  /* t1->Branch("Proton_Energy_RF_GeV",                     &fProton_Energy_RF_GeV,                     "fProton_Energy_RF_GeV/D"); */
  /* t1->Branch("Proton_Mom_RF_GeV",                        &fProton_Mom_RF_GeV,                        "fProton_Mom_RF_GeV/D"); */
  /* t1->Branch("Proton_MomZ_RF_GeV",                       &fProton_MomZ_RF_GeV,                       "fProton_MomZ_RF_GeV/D"); */
  /* t1->Branch("Proton_MomX_RF_GeV",                       &fProton_MomX_RF_GeV,                       "fProton_MomX_RF_GeV/D"); */
  /* t1->Branch("Proton_MomY_RF_GeV",                       &fProton_MomY_RF_GeV,                       "fProton_MomY_RF_GeV/D"); */

  /* t1->Branch("Electron_Theta_RF",                        &fElectron_Theta_RF,                        "fElectron_Theta_RF/D"); */
  /* t1->Branch("Electron_Phi_RF",                          &fElectron_Phi_RF,                          "fElectron_Phi_RF/D"); */
  /* t1->Branch("Electron_Energy_RF_GeV",                   &fElectron_Energy_RF_GeV,                   "fElectron_Energy_RF_GeV/D"); */
  /* t1->Branch("Electron_Mom_RF_GeV",                      &fElectron_Mom_RF_GeV,                      "fElectron_Mom_RF_GeV/D"); */
  /* t1->Branch("Electron_MomX_RF_GeV",                     &fElectron_MomX_RF_GeV,                     "fElectron_MomX_RF_GeV/D"); */
  /* t1->Branch("Electron_MomY_RF_GeV",                     &fElectron_MomY_RF_GeV,                     "fElectron_MomY_RF_GeV/D"); */
  /* t1->Branch("Electron_MomZ_RF_GeV",                     &fElectron_MomZ_RF_GeV,                     "fElectron_MomZ_RF_GeV/D"); */

  /* t1->Branch("ScatElec_Theta_RF",                        &fScatElec_Theta_RF,                        "fScatElec_Theta_RF/D"); */
  /* t1->Branch("ScatElec_Phi_RF",                          &fScatElec_Phi_RF,                          "fScatElec_Phi_RF/D"); */
  /* t1->Branch("ScatElec_Energy_RF_GeV",                   &fScatElec_Energy_RF_GeV,                   "fScatElec_Energy_RF_GeV/D"); */
  /* t1->Branch("ScatElec_Mom_RF_GeV",                      &fScatElec_Mom_RF_GeV,                      "fScatElec_Mom_RF_GeV/D"); */
  /* t1->Branch("ScatElec_MomX_RF_GeV",                     &fScatElec_MomX_RF_GeV,                     "fScatElec_MomX_RF_GeV/D"); */
  /* t1->Branch("ScatElec_MomY_RF_GeV",                     &fScatElec_MomY_RF_GeV,                     "fScatElec_MomY_RF_GeV/D"); */
  /* t1->Branch("ScatElec_MomZ_RF_GeV",                     &fScatElec_MomZ_RF_GeV,                     "fScatElec_MomZ_RF_GeV/D"); */

  /* t1->Branch("Photon_Energy_RF",                     &fPhoton_Energy_RF,                     "fPhoton_Energy_RF/D"); */
  /* t1->Branch("Photon_Mom_RF",                        &fPhoton_Mom_RF,                        "fPhoton_Mom_RF/D"); */
  /* t1->Branch("Photon_MomX_RF",                       &fPhoton_MomX_RF,                       "fPhoton_MomX_RF/D"); */
  /* t1->Branch("Photon_MomY_RF",                       &fPhoton_MomY_RF,                       "fPhoton_MomY_RF/D"); */
  /* t1->Branch("Photon_MomZ_RF",                       &fPhoton_MomZ_RF,                       "fPhoton_MomZ_RF/D"); */

  /* t1->Branch("Photon_Theta_RF",                          &fPhoton_Theta_RF,                          "fPhoton_Theta_RF/D"); */
  /* t1->Branch("Photon_Phi_RF",                            &fPhoton_Phi_RF,                            "fPhoton_Phi_RF/D"); */
  /* t1->Branch("Photon_Energy_RF_GeV",                     &fPhoton_Energy_RF_GeV,                     "fPhoton_Energy_RF_GeV/D"); */
  /* t1->Branch("Photon_Mom_RF_GeV",                        &fPhoton_Mom_RF_GeV,                        "fPhoton_Mom_RF_GeV/D"); */
  /* t1->Branch("Photon_MomX_RF_GeV",                       &fPhoton_MomX_RF_GeV,                       "fPhoton_MomX_RF_GeV/D"); */
  /* t1->Branch("Photon_MomY_RF_GeV",                       &fPhoton_MomY_RF_GeV,                       "fPhoton_MomY_RF_GeV/D"); */
  /* t1->Branch("Photon_MomZ_RF_GeV",                       &fPhoton_MomZ_RF_GeV,                       "fPhoton_MomZ_RF_GeV/D"); */

  /* t1->Branch("Pion_Theta_RF",                            &fPion_Theta_RF,                            "fPion_Theta_RF/D"); */
  /* t1->Branch("Pion_Phi_RF",                              &fPion_Phi_RF,                              "fPion_Phi_RF/D"); */
  /* t1->Branch("Pion_Energy_RF_GeV",                       &fPion_Energy_RF_GeV,                       "fPion_Energy_RF_GeV/D"); */
  /* t1->Branch("Pion_Mom_RF_GeV",                          &fPion_Mom_RF_GeV,                          "fPion_Mom_RF_GeV/D"); */
  /* t1->Branch("Pion_MomX_RF_GeV",                         &fPion_MomX_RF_GeV,                         "fPion_MomX_RF_GeV/D"); */
  /* t1->Branch("Pion_MomY_RF_GeV",                         &fPion_MomY_RF_GeV,                         "fPion_MomY_RF_GeV/D"); */
  /* t1->Branch("Pion_MomZ_RF_GeV",                         &fPion_MomZ_RF_GeV,                         "fPion_MomZ_RF_GeV/D"); */

  /* t1->Branch("Neutron_Theta_RF",                         &fNeutron_Theta_RF,                         "fNeutron_Theta_RF/D"); */
  /* t1->Branch("Neutron_Phi_RF",                           &fNeutron_Phi_RF,                           "fNeutron_Phi_RF/D"); */
  /* t1->Branch("Neutron_Energy_RF_GeV",                    &fNeutron_Energy_RF_GeV,                    "fNeutron_Energy_RF_GeV/D"); */
  /* t1->Branch("Neutron_Mom_RF_GeV",                       &fNeutron_Mom_RF_GeV,                       "fNeutron_Mom_RF_GeV/D"); */
  /* t1->Branch("Neutron_MomX_RF_GeV",                      &fNeutron_MomX_RF_GeV,                      "fNeutron_MomX_RF_GeV/D"); */
  /* t1->Branch("Neutron_MomY_RF_GeV",                      &fNeutron_MomY_RF_GeV,                      "fNeutron_MomY_RF_GeV/D"); */
  /* t1->Branch("Neutron_MomZ_RF_GeV",                      &fNeutron_MomZ_RF_GeV,                      "fNeutron_MomZ_RF_GeV/D"); */

  // ----------------------------------------------------------------------------------------------------------------------------------

  /* t1->Branch("Photon_Theta_RF",                           &fPhoton_Theta_RF,                           "fPhoton_Theta_RF/D"); */
  /* t1->Branch("Photon_Theta_Col",                          &fPhoton_Theta_Col,                          "fPhoton_Theta_Col/D"); */
  /* t1->Branch("NRecorded",                                 &fNRecorded,                                 "fNRecorded/I"); */
  /* t1->Branch("NGenerated",                                &fNGenerated,                                "fNGeneratedd/I"); */

  t1->Branch("Epsilon",                                   &fEpsilon,                                   "fEpsilon/D");
  t1->Branch("Phi",                                       &fPhi,                                       "fPhi/D");
  t1->Branch("PhiS",                                      &fPhiS,                                      "fPhiS/D");
  t1->Branch("W_GeV",                                     &fW_GeV,                                     "fW_GeV/D");
  t1->Branch("W_Prime_GeV",                               &fW_Prime_GeV,                               "fW_Prime_GeV/D");

  t1->Branch("Qsq_GeV",                                   &fQsq_GeV,                                   "fQsq_GeV/D");
  t1->Branch("T_Para_GeV",                                &fT_Para_GeV,                                "fT_Para_GeV/D");
  t1->Branch("T_GeV",                                     &fT_GeV,                                     "fT_GeV/D");
  t1->Branch("x",                                         &fx,                                         "fx/D");
  t1->Branch("y",                                         &fy,                                         "fy/D");
  t1->Branch("z",                                         &fz,                                         "fz/D");

  t1->Branch("Flux_Factor_RF",                            &fFlux_Factor_RF,                            "fFlux_Factor_RF/D");
  t1->Branch("Flux_Factor_Col",                           &fFlux_Factor_Col,                           "fFlux_Factor_Col/D");
  t1->Branch("Jacobian_CM",                               &fJacobian_CM,                               "fJacobian_CM/D");
  t1->Branch("Jacobian_CM_RF",                            &fJacobian_CM_RF,                            "fJacobian_CM_RF/D");
  t1->Branch("Jacobian_CM_Col",                           &fJacobian_CM_Col,                           "fJacobian_CM_Col/D");
  t1->Branch("EventWeight",                               &fEventWeight,                               "fEventWeight/D");
  t1->Branch("Sigma_Col",                                 &fSigma_Col,                                 "fSigma_Col/D");
  t1->Branch("Sig_VR"  ,                                  &fSig_VR,                                    "fSig_VR/D");
  t1->Branch("Sig_L",                                     &fSig_L,                                     "fSig_L/D");
  t1->Branch("Sig_T",                                     &fSig_T,                                     "fSig_T/D");

  t1->Branch("A",                                         &fA,                                         "fA/D");
  t1->Branch("Vertex_X",                                  &fVertex_X,                                  "fVertex_X/D");
  t1->Branch("Vertex_Y",                                  &fVertex_Y,                                  "fVertex_Y/D");
  t1->Branch("Vertex_Z",                                  &fVertex_Z,                                  "fVertex_Z/D");
  
  /* t1->Branch("Pion_Energy_CM_GeV",                        &fPion_Energy_CM_GeV,                        "fPion_Energy_CM_GeV/D"); */
  /* t1->Branch("Pion_Mom_CM_GeV",                           &fPion_Mom_CM_GeV,                           "fPion_Mom_CM_GeV/D"); */
  /* t1->Branch("BetaX_Col_RF",                              &fBetaX_Col_RF,                              "fBetaX_Col_RF/D"); */
  /* t1->Branch("BetaY_Col_RF",                              &fBetaY_Col_RF,                              "fBetaY_Col_RF/D"); */
  /* t1->Branch("BetaZ_Col_RF",                              &fBetaZ_Col_RF,                              "fBetaZ_Col_RF/D"); */
  /* t1->Branch("Beta_Col_RF",                               &fBeta_Col_RF,                               "fBeta_Col_RF/D"); */
  /* t1->Branch("Gamma_Col_RF",                              &fGamma_Col_RF,                              "fGamma_Col_RF/D"); */
  /* t1->Branch("Beta_CM_RF",                                &fBeta_CM_RF,                                "fBeta_CM_RF/D"); */
  /* t1->Branch("Gamma_CM_RF",                               &fGamma_CM_RF,                               "fGamma_CM_RF/D"); */

  /* t1->Branch("XMomConserve",                              &fXMomConserve,                              "fXMomConserve/D"); */
  /* t1->Branch("YMomConserve",                              &fYMomConserve,                              "fYMomConserve/D"); */
  /* t1->Branch("ZMomConserve",                              &fZMomConserve,                              "fZMomConserve/D"); */
  /* t1->Branch("EnergyConserve",                            &fEnergyConserve,                            "fEnergyConserve/D"); */

  /* t1->Branch("XMomConserve_RF",                           &fXMomConserve_RF,                           "fXMomConserve_RF/D"); */
  /* t1->Branch("YMomConserve_RF",                           &fYMomConserve_RF,                           "fYMomConserve_RF/D"); */
  /* t1->Branch("ZMomConserve_RF",                           &fZMomConserve_RF,                           "fZMomConserve_RF/D"); */
  /* t1->Branch("EnergyConserve_RF",                         &fEnergyConserve_RF,                         "fEnergyConserve_RF/D"); */
  /* t1->Branch("testsig",                                   &ftestsig,                                   "ftestsig/D"); */

}
