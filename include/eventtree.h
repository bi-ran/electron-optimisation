#ifndef EVENTTREE_H
#define EVENTTREE_H

#include "TTree.h"

#include <vector>

#define BRANCHES(ACTION)                              \
   BRANCHESMC(ACTION)                                 \
   BRANCHESDATA(ACTION)                               \

#define BRANCHESMC(ACTION)                            \
   ACTION(Int_t, nPUInfo)                             \
   ACTION(std::vector<int>*, nPU)                     \
   ACTION(std::vector<int>*, puBX)                    \
   ACTION(std::vector<float>*, puTrue)                \
   ACTION(Int_t, nMC)                                 \
   ACTION(std::vector<int>*, mcPID)                   \
   ACTION(std::vector<int>*, mcStatus)                \
   ACTION(std::vector<float>*, mcVtx_x)               \
   ACTION(std::vector<float>*, mcVtx_y)               \
   ACTION(std::vector<float>*, mcVtx_z)               \
   ACTION(std::vector<float>*, mcPt)                  \
   ACTION(std::vector<float>*, mcEta)                 \
   ACTION(std::vector<float>*, mcPhi)                 \
   ACTION(std::vector<float>*, mcE)                   \
   ACTION(std::vector<float>*, mcEt)                  \
   ACTION(std::vector<float>*, mcMass)                \
   ACTION(std::vector<int>*, mcParentage)             \
   ACTION(std::vector<int>*, mcMomPID)                \
   ACTION(std::vector<float>*, mcMomPt)               \
   ACTION(std::vector<float>*, mcMomEta)              \
   ACTION(std::vector<float>*, mcMomPhi)              \
   ACTION(std::vector<float>*, mcMomMass)             \
   ACTION(std::vector<int>*, mcGMomPID)               \
   ACTION(std::vector<int>*, mcIndex)                 \
   ACTION(std::vector<float>*, mcCalIsoDR03)          \
   ACTION(std::vector<float>*, mcCalIsoDR04)          \
   ACTION(std::vector<float>*, mcTrkIsoDR03)          \
   ACTION(std::vector<float>*, mcTrkIsoDR04)          \

#define BRANCHESDATA(ACTION)                          \
   ACTION(UInt_t, run)                                \
   ACTION(ULong64_t, event)                           \
   ACTION(UInt_t, lumis)                              \
   ACTION(Bool_t, isData)                             \
   ACTION(Int_t, nEle)                                \
   ACTION(std::vector<int>*, eleCharge)               \
   ACTION(std::vector<int>*, eleChargeConsistent)     \
   ACTION(std::vector<int>*, eleSCPixCharge)          \
   ACTION(std::vector<int>*, eleCtfCharge)            \
   ACTION(std::vector<float>*, eleEn)                 \
   ACTION(std::vector<float>*, eleD0)                 \
   ACTION(std::vector<float>*, eleDz)                 \
   ACTION(std::vector<float>*, eleD0Err)              \
   ACTION(std::vector<float>*, eleDzErr)              \
   ACTION(std::vector<float>*, eleTrkPt)              \
   ACTION(std::vector<float>*, eleTrkEta)             \
   ACTION(std::vector<float>*, eleTrkPhi)             \
   ACTION(std::vector<int>*, eleTrkCharge)            \
   ACTION(std::vector<float>*, eleTrkChi2)            \
   ACTION(std::vector<float>*, eleTrkNdof)            \
   ACTION(std::vector<float>*, eleTrkNormalizedChi2)  \
   ACTION(std::vector<int>*, eleTrkValidHits)         \
   ACTION(std::vector<int>*, eleTrkLayers)            \
   ACTION(std::vector<float>*, elePt)                 \
   ACTION(std::vector<float>*, eleEta)                \
   ACTION(std::vector<float>*, elePhi)                \
   ACTION(std::vector<float>*, eleSCEn)               \
   ACTION(std::vector<float>*, eleESEn)               \
   ACTION(std::vector<float>*, eleSCEta)              \
   ACTION(std::vector<float>*, eleSCPhi)              \
   ACTION(std::vector<float>*, eleSCRawEn)            \
   ACTION(std::vector<float>*, eleSCEtaWidth)         \
   ACTION(std::vector<float>*, eleSCPhiWidth)         \
   ACTION(std::vector<float>*, eleHoverE)             \
   ACTION(std::vector<float>*, eleHoverEBc)           \
   ACTION(std::vector<float>*, eleEoverP)             \
   ACTION(std::vector<float>*, eleEoverPInv)          \
   ACTION(std::vector<float>*, eleBrem)               \
   ACTION(std::vector<float>*, eledEtaAtVtx)          \
   ACTION(std::vector<float>*, eledPhiAtVtx)          \
   ACTION(std::vector<float>*, eleSigmaIEtaIEta)      \
   ACTION(std::vector<float>*, eleSigmaIEtaIEta_2012) \
   ACTION(std::vector<float>*, eleSigmaIPhiIPhi)      \
   ACTION(std::vector<int>*, eleMissHits)             \
   ACTION(std::vector<float>*, eleESEffSigmaRR)       \
   ACTION(std::vector<float>*, elePFChIso)            \
   ACTION(std::vector<float>*, elePFPhoIso)           \
   ACTION(std::vector<float>*, elePFNeuIso)           \
   ACTION(std::vector<float>*, elePFPUIso)            \
   ACTION(std::vector<float>*, elePFChIso03)          \
   ACTION(std::vector<float>*, elePFPhoIso03)         \
   ACTION(std::vector<float>*, elePFNeuIso03)         \
   ACTION(std::vector<float>*, elePFChIso04)          \
   ACTION(std::vector<float>*, elePFPhoIso04)         \
   ACTION(std::vector<float>*, elePFNeuIso04)         \
   ACTION(std::vector<float>*, eleR9)                 \
   ACTION(std::vector<float>*, eleE3x3)               \
   ACTION(std::vector<float>*, eleE5x5)               \
   ACTION(std::vector<float>*, eleR9Full5x5)          \
   ACTION(std::vector<float>*, eleE3x3Full5x5)        \
   ACTION(std::vector<float>*, eleE5x5Full5x5)        \
   ACTION(std::vector<int>*, NClusters)               \
   ACTION(std::vector<int>*, NEcalClusters)           \
   ACTION(std::vector<float>*, eleSeedEn)             \
   ACTION(std::vector<float>*, eleSeedEta)            \
   ACTION(std::vector<float>*, eleSeedPhi)            \
   ACTION(std::vector<float>*, eleSeedCryEta)         \
   ACTION(std::vector<float>*, eleSeedCryPhi)         \
   ACTION(std::vector<float>*, eleSeedCryIeta)        \
   ACTION(std::vector<float>*, eleSeedCryIphi)        \
   ACTION(std::vector<float>*, eleBC1E)               \
   ACTION(std::vector<float>*, eleBC1Eta)             \
   ACTION(std::vector<float>*, eleBC2E)               \
   ACTION(std::vector<float>*, eleBC2Eta)             \
   ACTION(std::vector<int>*, eleIDVeto)               \
   ACTION(std::vector<int>*, eleIDLoose)              \
   ACTION(std::vector<int>*, eleIDMedium)             \
   ACTION(std::vector<int>*, eleIDTight)              \
   ACTION(std::vector<int>*, elepassConversionVeto)   \
   ACTION(std::vector<float>*, eleEffAreaTimesRho)    \

#define ZERO(type, var) var = 0;
#define DECLARE(type, var) type var;
#define READ(type, var)                               \
   t->SetBranchStatus(#var, 1);                       \
   t->SetBranchAddress(#var, &var);                   \

class eventtree {
   public:
      eventtree() { BRANCHES(ZERO) };
      eventtree(TTree* t, bool isdata) : eventtree() { read(t, isdata); }
      ~eventtree() { };

      void read(TTree* t, bool isdata) {
         if (!isdata) {
            BRANCHESMC(READ) }
         BRANCHESDATA(READ)
      };

      BRANCHES(DECLARE)
};

#undef BRANCHES
#undef ZERO
#undef DECLARE
#undef READ

#endif  /* EVENTTREE_H */
