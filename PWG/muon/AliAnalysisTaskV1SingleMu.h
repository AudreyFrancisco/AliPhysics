#ifndef ALIANALYSISTASKV1SINGLEMU_H
#define ALIANALYSISTASKV1SINGLEMU_H

//
// AliAnalysisTaskV1SingleMu
// Analysis task for v1 of single muons in the spectrometer with the scalar product method
//
//  Author: Audrey Francisco
//
//
//

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisMuonUtility.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
// #include "AliMuonTrackCuts.h"

class TObjArray;
class THnSparse;
class AliMergeableCollection;
class AliMuonTrackCuts;
class AliMuonEventCuts;
class AliUtilityMuonAncestor;


class AliAnalysisTaskV1SingleMu : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskV1SingleMu();
  AliAnalysisTaskV1SingleMu(const char *name, const AliMuonTrackCuts& cuts);
  virtual ~AliAnalysisTaskV1SingleMu();

  /// Get muon track cuts
  AliMuonTrackCuts* GetMuonTrackCuts() { return fMuonTrackCuts; }
  AliMuonEventCuts* GetMuonEventCuts() { return fMuonEventCuts; }

  void UserCreateOutputObjects();
  void UserExec(Option_t *option);
  // void ProcessEvent(TString physSel, const TObjArray& selectTrigClasses, TString centrality);
  void  Terminate(Option_t *option);

  void NotifyRun();

 private:
  TObject* GetMergeableObject ( TString identifier, TString objectName );

  AliAnalysisTaskV1SingleMu(const AliAnalysisTaskV1SingleMu&);
  AliAnalysisTaskV1SingleMu& operator=(const AliAnalysisTaskV1SingleMu&);

  Int_t GetParticleType ( AliVParticle* track );
  
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

  AliMuonEventCuts* fMuonEventCuts; ///< Muon event cuts
  AliMuonTrackCuts* fMuonTrackCuts;  ///< Muon event cuts
  AliESDEvent* fESDEvent;      //!< ESD event, not owner
  AliAODEvent* fAODEvent;      //!< AOD event, not owner
  AliUtilityMuonAncestor* fUtilityMuonAncestor; ///< Utility to get the muon ancestor for MC
  Bool_t fCutOnDimu;
  Int_t fNPtBins;
  Int_t fHarmonic;
  TString fNormMethod;
  AliMergeableCollection* fMergeableCollection; //!<! collection of mergeable objects
  THnSparse* fSparse; ///< CF container

  ClassDef(AliAnalysisTaskV1SingleMu, 1); // Single muon analysis
};

#endif
