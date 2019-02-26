#ifndef ELECTRONTREE_H
#define ELECTRONTREE_H

#include "TTree.h"

#include <vector>
#include <algorithm>
#include <iterator>

#include "eventtree.h"
#include "l1tree.h"

#define BRANCHES(ACTION)                              \
   VARBRANCHES(ACTION)                                \
   NEWVARBRANCHES(ACTION)                             \
   VECBRANCHESMC(ACTION)                              \
   VECBRANCHESDATA(ACTION)                            \
   NEWVECBRANCHES(ACTION)                             \
   L1VARBRANCHES(ACTION)                              \
   L1VECBRANCHES(ACTION)                              \

#define VARBRANCHES(ACTION)                           \
   ACTION(UInt_t, run)                                \
   ACTION(ULong64_t, event)                           \
   ACTION(UInt_t, lumis)                              \
   ACTION(Bool_t, isData)                             \
   ACTION(Int_t, nPUInfo)                             \
   ACTION(Int_t, nMC)                                 \
   ACTION(Int_t, nEle)                                \

#define NEWVARBRANCHES(ACTION)                        \
   ACTION(int, hiBin)                                 \
   ACTION(float, hiHF)                                \
   ACTION(float, ncoll)                               \

#define VECBRANCHESMC(ACTION)                         \
   ACTION(std::vector<int>, nPU)                      \
   ACTION(std::vector<int>, puBX)                     \
   ACTION(std::vector<float>, puTrue)                 \
   ACTION(std::vector<int>, mcPID)                    \
   ACTION(std::vector<int>, mcStatus)                 \
   ACTION(std::vector<float>, mcVtx_x)                \
   ACTION(std::vector<float>, mcVtx_y)                \
   ACTION(std::vector<float>, mcVtx_z)                \
   ACTION(std::vector<float>, mcPt)                   \
   ACTION(std::vector<float>, mcEta)                  \
   ACTION(std::vector<float>, mcPhi)                  \
   ACTION(std::vector<float>, mcE)                    \
   ACTION(std::vector<float>, mcEt)                   \
   ACTION(std::vector<float>, mcMass)                 \
   ACTION(std::vector<int>, mcParentage)              \
   ACTION(std::vector<int>, mcMomPID)                 \
   ACTION(std::vector<float>, mcMomPt)                \
   ACTION(std::vector<float>, mcMomEta)               \
   ACTION(std::vector<float>, mcMomPhi)               \
   ACTION(std::vector<float>, mcMomMass)              \
   ACTION(std::vector<int>, mcGMomPID)                \
   ACTION(std::vector<int>, mcIndex)                  \
   ACTION(std::vector<float>, mcCalIsoDR03)           \
   ACTION(std::vector<float>, mcCalIsoDR04)           \
   ACTION(std::vector<float>, mcTrkIsoDR03)           \
   ACTION(std::vector<float>, mcTrkIsoDR04)           \

#define VECBRANCHESDATA(ACTION)                       \
   ACTION(std::vector<int>, eleCharge)                \
   ACTION(std::vector<int>, eleChargeConsistent)      \
   ACTION(std::vector<int>, eleSCPixCharge)           \
   ACTION(std::vector<int>, eleCtfCharge)             \
   ACTION(std::vector<float>, eleEn)                  \
   ACTION(std::vector<float>, eleD0)                  \
   ACTION(std::vector<float>, eleDz)                  \
   ACTION(std::vector<float>, eleD0Err)               \
   ACTION(std::vector<float>, eleDzErr)               \
   ACTION(std::vector<float>, eleTrkPt)               \
   ACTION(std::vector<float>, eleTrkEta)              \
   ACTION(std::vector<float>, eleTrkPhi)              \
   ACTION(std::vector<int>, eleTrkCharge)             \
   ACTION(std::vector<float>, eleTrkChi2)             \
   ACTION(std::vector<float>, eleTrkNdof)             \
   ACTION(std::vector<float>, eleTrkNormalizedChi2)   \
   ACTION(std::vector<int>, eleTrkValidHits)          \
   ACTION(std::vector<int>, eleTrkLayers)             \
   ACTION(std::vector<float>, elePt)                  \
   ACTION(std::vector<float>, eleEta)                 \
   ACTION(std::vector<float>, elePhi)                 \
   ACTION(std::vector<float>, eleSCEn)                \
   ACTION(std::vector<float>, eleESEn)                \
   ACTION(std::vector<float>, eleSCEta)               \
   ACTION(std::vector<float>, eleSCPhi)               \
   ACTION(std::vector<float>, eleSCRawEn)             \
   ACTION(std::vector<float>, eleSCEtaWidth)          \
   ACTION(std::vector<float>, eleSCPhiWidth)          \
   ACTION(std::vector<float>, eleHoverE)              \
   ACTION(std::vector<float>, eleHoverEBc)            \
   ACTION(std::vector<float>, eleEoverP)              \
   ACTION(std::vector<float>, eleEoverPInv)           \
   ACTION(std::vector<float>, eleBrem)                \
   ACTION(std::vector<float>, eledEtaAtVtx)           \
   ACTION(std::vector<float>, eledPhiAtVtx)           \
   ACTION(std::vector<float>, eleSigmaIEtaIEta)       \
   ACTION(std::vector<float>, eleSigmaIEtaIEta_2012)  \
   ACTION(std::vector<float>, eleSigmaIPhiIPhi)       \
   ACTION(std::vector<int>, eleMissHits)              \
   ACTION(std::vector<float>, eleESEffSigmaRR)        \
   ACTION(std::vector<float>, elePFChIso)             \
   ACTION(std::vector<float>, elePFPhoIso)            \
   ACTION(std::vector<float>, elePFNeuIso)            \
   ACTION(std::vector<float>, elePFPUIso)             \
   ACTION(std::vector<float>, elePFChIso03)           \
   ACTION(std::vector<float>, elePFPhoIso03)          \
   ACTION(std::vector<float>, elePFNeuIso03)          \
   ACTION(std::vector<float>, elePFChIso04)           \
   ACTION(std::vector<float>, elePFPhoIso04)          \
   ACTION(std::vector<float>, elePFNeuIso04)          \
   ACTION(std::vector<float>, eleR9)                  \
   ACTION(std::vector<float>, eleE3x3)                \
   ACTION(std::vector<float>, eleE5x5)                \
   ACTION(std::vector<float>, eleR9Full5x5)           \
   ACTION(std::vector<float>, eleE3x3Full5x5)         \
   ACTION(std::vector<float>, eleE5x5Full5x5)         \
   ACTION(std::vector<int>, NClusters)                \
   ACTION(std::vector<int>, NEcalClusters)            \
   ACTION(std::vector<float>, eleSeedEn)              \
   ACTION(std::vector<float>, eleSeedEta)             \
   ACTION(std::vector<float>, eleSeedPhi)             \
   ACTION(std::vector<float>, eleSeedCryEta)          \
   ACTION(std::vector<float>, eleSeedCryPhi)          \
   ACTION(std::vector<float>, eleSeedCryIeta)         \
   ACTION(std::vector<float>, eleSeedCryIphi)         \
   ACTION(std::vector<float>, eleBC1E)                \
   ACTION(std::vector<float>, eleBC1Eta)              \
   ACTION(std::vector<float>, eleBC2E)                \
   ACTION(std::vector<float>, eleBC2Eta)              \
   ACTION(std::vector<int>, eleIDVeto)                \
   ACTION(std::vector<int>, eleIDLoose)               \
   ACTION(std::vector<int>, eleIDMedium)              \
   ACTION(std::vector<int>, eleIDTight)               \
   ACTION(std::vector<int>, elepassConversionVeto)    \
   ACTION(std::vector<float>, eleEffAreaTimesRho)     \

#define NEWVECBRANCHES(ACTION)                        \
   ACTION(std::vector<int>, eleGenMatchIndex)         \
   ACTION(std::vector<int>, mcRecoMatchIndex)         \
   ACTION(std::vector<int>, hlt)                      \
   ACTION(std::vector<double>, pt_e20)                \
   ACTION(std::vector<double>, eta_e20)               \
   ACTION(std::vector<double>, phi_e20)               \
   ACTION(std::vector<int>, n_e20)                    \
   ACTION(std::vector<double>, pt_e10e10m50)          \
   ACTION(std::vector<double>, eta_e10e10m50)         \
   ACTION(std::vector<double>, phi_e10e10m50)         \
   ACTION(std::vector<int>, n_e10e10m50)              \

#define L1VARBRANCHES(ACTION)                         \
   ACTION(short, nEGs)                                \

#define L1VECBRANCHES(ACTION)                         \
   ACTION(std::vector<float>, egEt)                   \
   ACTION(std::vector<float>, egEta)                  \
   ACTION(std::vector<float>, egPhi)                  \
   ACTION(std::vector<short>, egIEt)                  \
   ACTION(std::vector<short>, egIEta)                 \
   ACTION(std::vector<short>, egIPhi)                 \
   ACTION(std::vector<short>, egIso)                  \
   ACTION(std::vector<short>, egBx)                   \
   ACTION(std::vector<short>, egTowerIPhi)            \
   ACTION(std::vector<short>, egTowerIEta)            \
   ACTION(std::vector<short>, egRawEt)                \
   ACTION(std::vector<short>, egIsoEt)                \
   ACTION(std::vector<short>, egFootprintEt)          \
   ACTION(std::vector<short>, egNTT)                  \
   ACTION(std::vector<short>, egShape)                \
   ACTION(std::vector<short>, egTowerHoE)             \
   ACTION(std::vector<short>, egHwQual)               \

#define INVALID(type, var) var = -1;
#define DECLARE(type, var) type var;
#define CREATE(type, var) t->Branch(#var, &var);
#define CLEAR(type, var) var.clear();
#define VARCOPY(type, var) var = evtt->var;
#define VECCOPY(type, var)                                                    \
   std::copy(evtt->var->begin(), evtt->var->end(), std::back_inserter(var));
#define L1VARCOPY(type, var) var = l1ot->var;
#define L1VECCOPY(type, var)                                                  \
   std::copy(l1ot->var->begin(), l1ot->var->end(), std::back_inserter(var));

class electrontree {
   public:
      electrontree() { VARBRANCHES(INVALID) };
      electrontree(TTree* t, bool isdata)
         : electrontree() { this->isdata = isdata; branch(t); };
      ~electrontree() {};

      void branch(TTree* t) {
         VARBRANCHES(CREATE)
         if (!isdata) {
            VECBRANCHESMC(CREATE) }
         else {
            L1VARBRANCHES(CREATE)
            L1VECBRANCHES(CREATE) }
         VECBRANCHESDATA(CREATE)
         NEWVARBRANCHES(CREATE)
         NEWVECBRANCHES(CREATE)
      };

      void clear() {
         VECBRANCHESMC(CLEAR)
         VECBRANCHESDATA(CLEAR)
         NEWVECBRANCHES(CLEAR)
         L1VECBRANCHES(CLEAR)
      };

      void copy(eventtree* evtt) {
         VARBRANCHES(VARCOPY)
         if (!isdata) {
            VECBRANCHESMC(VECCOPY) }
         VECBRANCHESDATA(VECCOPY)
      };

      void copy(l1tree* l1ot) {
         if (isdata) {
            L1VARBRANCHES(L1VARCOPY)
            L1VECBRANCHES(L1VECCOPY)
         }
      };

      BRANCHES(DECLARE)

   private:
      bool isdata;
};

#endif  /* ELECTRONTREE_H */
