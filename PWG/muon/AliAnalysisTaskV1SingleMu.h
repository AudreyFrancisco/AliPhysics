#ifndef ALIANALYSISTASKV1SINGLEMU_H
#define ALIANALYSISTASKV1SINGLEMU_H

//
// AliAnalysisTaskV1SingleMu
// Analysis task for v1 of single muons in the spectrometer with the scalar product method
//
//  Author: Audrey Francisco from AliAnalysisDimu
//


#include "AliAnalysisTaskSE.h"
#include "AliAnalysisMuonUtility.h"
// #include "AliAODEvent.h"
// #include "AliESDEvent.h"
#include "AliMuonTrackCuts.h"
#include "AliMuonEventCuts.h"
// #include "AliMuonTrackCuts.h"

class TObjArray;
class THnSparse;
class AliMergeableCollection;
class AliUtilityMuonAncestor;


class AliAnalysisTaskV1SingleMu : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskV1SingleMu();
  AliAnalysisTaskV1SingleMu(const char *name);//, const AliMuonTrackCuts& cuts);
  virtual ~AliAnalysisTaskV1SingleMu();

  void UserCreateOutputObjects();
  void UserExec(Option_t *option);
  void NotifyRun();
  void  Terminate(Option_t *option);

  /// Get muon event cuts
  AliMuonEventCuts* GetMuonEventCuts() { return &fMuonEventCuts; }
  /// Get muon track cuts
  AliMuonTrackCuts* GetMuonTrackCuts() { return &fMuonTrackCuts; }

  /// Set muon event cuts
  void SetMuonEventCuts ( AliMuonEventCuts* muonEventCuts ) { fMuonEventCuts = *muonEventCuts; }
  /// Set muon track cuts
  void SetMuonTrackCuts ( AliMuonTrackCuts* muonTrackCuts ) { fMuonTrackCuts = *muonTrackCuts; }

  enum {
    kStepReconstructed,  ///< Reconstructed tracks
    kStepGeneratedMC,    ///< Generated tracks (MC)
    kNsteps              ///< Number of steps
  };
  enum {
    kHvarPt,         ///< Pt at vertex
    kHvarEta,        ///< Pseudo-Rapidity
    kHvarPhi,        ///< Phi
    kHvarCharge,     ///< Particle charge
    kHCentrality,    ///< event centrality
    //TODO :add SP var
    kNvars           ///< THnSparse dimensions
  };

 private:
  TObject* GetMergeableObject ( TString identifier, TString objectName );

  AliAnalysisTaskV1SingleMu(const AliAnalysisTaskV1SingleMu&);
  AliAnalysisTaskV1SingleMu& operator=(const AliAnalysisTaskV1SingleMu&);

  Int_t GetParticleType ( AliVParticle* track );

  //add enum for MC (AliUtilityMuonAncestor)
    enum {
    kCharmMu,       ///< Mu from charm
    kBeautyMu,      ///< Mu from beauty
    kQuarkoniumMu,  ///< Mu from resonance
    kWbosonMu,      ///< Mu from W
    kZbosonMu,      ///< Mu from Z
    kDecayMu,       ///< Decay mu
    kSecondaryMu,   ///< Secondary mu
    kRecoHadron,    ///< Reconstructed hadron
    kUnidentified,  ///< Particle that fails matching kine
    kNtrackSources  ///< Total number of track sources
  };

  AliMuonEventCuts fMuonEventCuts;  ///< Muon event cuts
  AliMuonTrackCuts fMuonTrackCuts;  ///< Muon event cuts
  // AliESDEvent* fESDEvent;      //!< ESD event, not owner
  // AliAODEvent* fAODEvent;      //!< AOD event, not owner
  AliUtilityMuonAncestor* fUtilityMuonAncestor; ///< Utility to get the muon ancestor for MC
  Int_t fNPtBins;
  Int_t fHarmonic;
  TString fNormMethod;
  AliMergeableCollection* fMergeableCollection; //!<! collection of mergeable objects
  THnSparse* fSparse; ///< CF container

  ClassDef(AliAnalysisTaskV1SingleMu, 1); // Single muon v1 analysis
};

#endif
