#ifndef EVENTTREE_H
#define EVENTTREE_H

#include "TTree.h"

#include <vector>

#include "common.h"

#define BRANCHES_EVT(ACTION, ...)                                       \
   B_ALE_M(ACTION, ## __VA_ARGS__)                                      \
   B_ALE_D(ACTION, ## __VA_ARGS__)                                      \

#define B_ALE_M(ACTION, ...)                                            \
   ACTION(Int_t, nMC, ## __VA_ARGS__)                                   \
   ACTION(std::vector<int>*, mcPID, ## __VA_ARGS__)                     \
   ACTION(std::vector<int>*, mcStatus, ## __VA_ARGS__)                  \
   ACTION(std::vector<float>*, mcVtx_x, ## __VA_ARGS__)                 \
   ACTION(std::vector<float>*, mcVtx_y, ## __VA_ARGS__)                 \
   ACTION(std::vector<float>*, mcVtx_z, ## __VA_ARGS__)                 \
   ACTION(std::vector<float>*, mcPt, ## __VA_ARGS__)                    \
   ACTION(std::vector<float>*, mcEta, ## __VA_ARGS__)                   \
   ACTION(std::vector<float>*, mcPhi, ## __VA_ARGS__)                   \
   ACTION(std::vector<float>*, mcE, ## __VA_ARGS__)                     \
   ACTION(std::vector<float>*, mcEt, ## __VA_ARGS__)                    \
   ACTION(std::vector<float>*, mcMass, ## __VA_ARGS__)                  \
   ACTION(std::vector<int>*, mcParentage, ## __VA_ARGS__)               \
   ACTION(std::vector<int>*, mcMomPID, ## __VA_ARGS__)                  \
   ACTION(std::vector<float>*, mcMomPt, ## __VA_ARGS__)                 \
   ACTION(std::vector<float>*, mcMomEta, ## __VA_ARGS__)                \
   ACTION(std::vector<float>*, mcMomPhi, ## __VA_ARGS__)                \
   ACTION(std::vector<float>*, mcMomMass, ## __VA_ARGS__)               \
   ACTION(std::vector<int>*, mcGMomPID, ## __VA_ARGS__)                 \
   ACTION(std::vector<int>*, mcIndex, ## __VA_ARGS__)                   \
   ACTION(std::vector<float>*, mcCalIsoDR03, ## __VA_ARGS__)            \
   ACTION(std::vector<float>*, mcCalIsoDR04, ## __VA_ARGS__)            \
   ACTION(std::vector<float>*, mcTrkIsoDR03, ## __VA_ARGS__)            \
   ACTION(std::vector<float>*, mcTrkIsoDR04, ## __VA_ARGS__)            \

#define B_ALE_D(ACTION, ...)                                            \
   ACTION(UInt_t, run, ## __VA_ARGS__)                                  \
   ACTION(ULong64_t, event, ## __VA_ARGS__)                             \
   ACTION(UInt_t, lumis, ## __VA_ARGS__)                                \
   ACTION(Int_t, nEle, ## __VA_ARGS__)                                  \
   ACTION(std::vector<int>*, eleCharge, ## __VA_ARGS__)                 \
   ACTION(std::vector<int>*, eleChargeConsistent, ## __VA_ARGS__)       \
   ACTION(std::vector<int>*, eleSCPixCharge, ## __VA_ARGS__)            \
   ACTION(std::vector<int>*, eleCtfCharge, ## __VA_ARGS__)              \
   ACTION(std::vector<float>*, eleEn, ## __VA_ARGS__)                   \
   ACTION(std::vector<float>*, eleD0, ## __VA_ARGS__)                   \
   ACTION(std::vector<float>*, eleDz, ## __VA_ARGS__)                   \
   ACTION(std::vector<float>*, eleD0Err, ## __VA_ARGS__)                \
   ACTION(std::vector<float>*, eleDzErr, ## __VA_ARGS__)                \
   ACTION(std::vector<float>*, eleTrkPt, ## __VA_ARGS__)                \
   ACTION(std::vector<float>*, eleTrkEta, ## __VA_ARGS__)               \
   ACTION(std::vector<float>*, eleTrkPhi, ## __VA_ARGS__)               \
   ACTION(std::vector<int>*, eleTrkCharge, ## __VA_ARGS__)              \
   ACTION(std::vector<float>*, eleTrkPtErr, ## __VA_ARGS__)             \
   ACTION(std::vector<float>*, eleTrkChi2, ## __VA_ARGS__)              \
   ACTION(std::vector<float>*, eleTrkNdof, ## __VA_ARGS__)              \
   ACTION(std::vector<float>*, eleTrkNormalizedChi2, ## __VA_ARGS__)    \
   ACTION(std::vector<int>*, eleTrkValidHits, ## __VA_ARGS__)           \
   ACTION(std::vector<int>*, eleTrkLayers, ## __VA_ARGS__)              \
   ACTION(std::vector<float>*, elePt, ## __VA_ARGS__)                   \
   ACTION(std::vector<float>*, eleEta, ## __VA_ARGS__)                  \
   ACTION(std::vector<float>*, elePhi, ## __VA_ARGS__)                  \
   ACTION(std::vector<float>*, eleSCEn, ## __VA_ARGS__)                 \
   ACTION(std::vector<float>*, eleESEn, ## __VA_ARGS__)                 \
   ACTION(std::vector<float>*, eleSCEta, ## __VA_ARGS__)                \
   ACTION(std::vector<float>*, eleSCPhi, ## __VA_ARGS__)                \
   ACTION(std::vector<float>*, eleSCRawEn, ## __VA_ARGS__)              \
   ACTION(std::vector<float>*, eleSCEtaWidth, ## __VA_ARGS__)           \
   ACTION(std::vector<float>*, eleSCPhiWidth, ## __VA_ARGS__)           \
   ACTION(std::vector<float>*, eleHoverE, ## __VA_ARGS__)               \
   ACTION(std::vector<float>*, eleHoverEBc, ## __VA_ARGS__)             \
   ACTION(std::vector<float>*, eleEoverP, ## __VA_ARGS__)               \
   ACTION(std::vector<float>*, eleEoverPInv, ## __VA_ARGS__)            \
   ACTION(std::vector<float>*, eleEcalE, ## __VA_ARGS__)                \
   ACTION(std::vector<float>*, elePAtVtx, ## __VA_ARGS__)               \
   ACTION(std::vector<float>*, elePAtSC, ## __VA_ARGS__)                \
   ACTION(std::vector<float>*, elePAtCluster, ## __VA_ARGS__)           \
   ACTION(std::vector<float>*, elePAtSeed, ## __VA_ARGS__)              \
   ACTION(std::vector<float>*, eleBrem, ## __VA_ARGS__)                 \
   ACTION(std::vector<float>*, eledEtaAtVtx, ## __VA_ARGS__)            \
   ACTION(std::vector<float>*, eledPhiAtVtx, ## __VA_ARGS__)            \
   ACTION(std::vector<float>*, eleSigmaIEtaIEta, ## __VA_ARGS__)        \
   ACTION(std::vector<float>*, eleSigmaIEtaIEta_2012, ## __VA_ARGS__)   \
   ACTION(std::vector<float>*, eleSigmaIPhiIPhi, ## __VA_ARGS__)        \
   ACTION(std::vector<int>*, eleMissHits, ## __VA_ARGS__)               \
   ACTION(std::vector<float>*, eleESEffSigmaRR, ## __VA_ARGS__)         \
   ACTION(std::vector<float>*, elePFChIso, ## __VA_ARGS__)              \
   ACTION(std::vector<float>*, elePFPhoIso, ## __VA_ARGS__)             \
   ACTION(std::vector<float>*, elePFNeuIso, ## __VA_ARGS__)             \
   ACTION(std::vector<float>*, elePFPUIso, ## __VA_ARGS__)              \
   ACTION(std::vector<float>*, elePFChIso03, ## __VA_ARGS__)            \
   ACTION(std::vector<float>*, elePFPhoIso03, ## __VA_ARGS__)           \
   ACTION(std::vector<float>*, elePFNeuIso03, ## __VA_ARGS__)           \
   ACTION(std::vector<float>*, elePFChIso04, ## __VA_ARGS__)            \
   ACTION(std::vector<float>*, elePFPhoIso04, ## __VA_ARGS__)           \
   ACTION(std::vector<float>*, elePFNeuIso04, ## __VA_ARGS__)           \
   ACTION(std::vector<float>*, eleR9, ## __VA_ARGS__)                   \
   ACTION(std::vector<float>*, eleE3x3, ## __VA_ARGS__)                 \
   ACTION(std::vector<float>*, eleE5x5, ## __VA_ARGS__)                 \
   ACTION(std::vector<float>*, eleR9Full5x5, ## __VA_ARGS__)            \
   ACTION(std::vector<float>*, eleE3x3Full5x5, ## __VA_ARGS__)          \
   ACTION(std::vector<float>*, eleE5x5Full5x5, ## __VA_ARGS__)          \
   ACTION(std::vector<int>*, NClusters, ## __VA_ARGS__)                 \
   ACTION(std::vector<int>*, NEcalClusters, ## __VA_ARGS__)             \
   ACTION(std::vector<float>*, eleSeedEn, ## __VA_ARGS__)               \
   ACTION(std::vector<float>*, eleSeedEta, ## __VA_ARGS__)              \
   ACTION(std::vector<float>*, eleSeedPhi, ## __VA_ARGS__)              \
   ACTION(std::vector<float>*, eleSeedCryEta, ## __VA_ARGS__)           \
   ACTION(std::vector<float>*, eleSeedCryPhi, ## __VA_ARGS__)           \
   ACTION(std::vector<float>*, eleSeedCryIeta, ## __VA_ARGS__)          \
   ACTION(std::vector<float>*, eleSeedCryIphi, ## __VA_ARGS__)          \
   ACTION(std::vector<float>*, eleBC1E, ## __VA_ARGS__)                 \
   ACTION(std::vector<float>*, eleBC1Eta, ## __VA_ARGS__)               \
   ACTION(std::vector<float>*, eleBC2E, ## __VA_ARGS__)                 \
   ACTION(std::vector<float>*, eleBC2Eta, ## __VA_ARGS__)               \
   ACTION(std::vector<int>*, eleIDVeto, ## __VA_ARGS__)                 \
   ACTION(std::vector<int>*, eleIDLoose, ## __VA_ARGS__)                \
   ACTION(std::vector<int>*, eleIDMedium, ## __VA_ARGS__)               \
   ACTION(std::vector<int>*, eleIDTight, ## __VA_ARGS__)                \
   ACTION(std::vector<float>*, eleEffAreaTimesRho, ## __VA_ARGS__)      \

class eventtree {
   public:
      eventtree() { do_mc_branches = 0; BRANCHES_EVT(ZERO) };
      eventtree(TTree* t, bool do_mc_branches)
         : eventtree() {
            this->do_mc_branches = do_mc_branches;
            read(t); }
      ~eventtree() { };

      void read(TTree* t) {
         if (do_mc_branches) {
            B_ALE_M(READ, t) }
         B_ALE_D(READ, t)
      };

      bool do_mc_branches;

      BRANCHES_EVT(DECLARE)
};

#endif  /* EVENTTREE_H */
