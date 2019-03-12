#ifndef ELECTRONTREE_H
#define ELECTRONTREE_H

#include "TTree.h"

#include <vector>
#include <algorithm>
#include <iterator>

#include "common.h"
#include "eventtree.h"
#include "l1tree.h"

#define B_VAR_D(ACTION, ...)                                            \
   ACTION(UInt_t, run, ## __VA_ARGS__)                                  \
   ACTION(ULong64_t, event, ## __VA_ARGS__)                             \
   ACTION(UInt_t, lumis, ## __VA_ARGS__)                                \
   ACTION(Int_t, nEle, ## __VA_ARGS__)                                  \

#define B_VAR_L(ACTION, ...)                                            \
   ACTION(short, nEGs, ## __VA_ARGS__)                                  \

#define B_VAR_M(ACTION, ...)                                            \
   ACTION(Int_t, nMC, ## __VA_ARGS__)                                   \

#define B_VAR_N(ACTION, ...)                                            \
   ACTION(int, hiBin, ## __VA_ARGS__)                                   \
   ACTION(float, hiHF, ## __VA_ARGS__)                                  \
   ACTION(float, ncoll, ## __VA_ARGS__)                                 \

#define B_VEC_D(ACTION, ...)                                            \
   ACTION(std::vector<int>, eleCharge, ## __VA_ARGS__)                  \
   ACTION(std::vector<int>, eleChargeConsistent, ## __VA_ARGS__)        \
   ACTION(std::vector<int>, eleSCPixCharge, ## __VA_ARGS__)             \
   ACTION(std::vector<int>, eleCtfCharge, ## __VA_ARGS__)               \
   ACTION(std::vector<float>, eleEn, ## __VA_ARGS__)                    \
   ACTION(std::vector<float>, eleD0, ## __VA_ARGS__)                    \
   ACTION(std::vector<float>, eleDz, ## __VA_ARGS__)                    \
   ACTION(std::vector<float>, eleD0Err, ## __VA_ARGS__)                 \
   ACTION(std::vector<float>, eleDzErr, ## __VA_ARGS__)                 \
   ACTION(std::vector<float>, eleTrkPt, ## __VA_ARGS__)                 \
   ACTION(std::vector<float>, eleTrkEta, ## __VA_ARGS__)                \
   ACTION(std::vector<float>, eleTrkPhi, ## __VA_ARGS__)                \
   ACTION(std::vector<int>, eleTrkCharge, ## __VA_ARGS__)               \
   ACTION(std::vector<float>, eleTrkPtErr, ## __VA_ARGS__)              \
   ACTION(std::vector<float>, eleTrkChi2, ## __VA_ARGS__)               \
   ACTION(std::vector<float>, eleTrkNdof, ## __VA_ARGS__)               \
   ACTION(std::vector<float>, eleTrkNormalizedChi2, ## __VA_ARGS__)     \
   ACTION(std::vector<int>, eleTrkValidHits, ## __VA_ARGS__)            \
   ACTION(std::vector<int>, eleTrkLayers, ## __VA_ARGS__)               \
   ACTION(std::vector<float>, elePt, ## __VA_ARGS__)                    \
   ACTION(std::vector<float>, eleEta, ## __VA_ARGS__)                   \
   ACTION(std::vector<float>, elePhi, ## __VA_ARGS__)                   \
   ACTION(std::vector<float>, eleSCEn, ## __VA_ARGS__)                  \
   ACTION(std::vector<float>, eleESEn, ## __VA_ARGS__)                  \
   ACTION(std::vector<float>, eleSCEta, ## __VA_ARGS__)                 \
   ACTION(std::vector<float>, eleSCPhi, ## __VA_ARGS__)                 \
   ACTION(std::vector<float>, eleSCRawEn, ## __VA_ARGS__)               \
   ACTION(std::vector<float>, eleSCEtaWidth, ## __VA_ARGS__)            \
   ACTION(std::vector<float>, eleSCPhiWidth, ## __VA_ARGS__)            \
   ACTION(std::vector<float>, eleHoverE, ## __VA_ARGS__)                \
   ACTION(std::vector<float>, eleHoverEBc, ## __VA_ARGS__)              \
   ACTION(std::vector<float>, eleEoverP, ## __VA_ARGS__)                \
   ACTION(std::vector<float>, eleEoverPInv, ## __VA_ARGS__)             \
   ACTION(std::vector<float>, eleEcalE, ## __VA_ARGS__)                 \
   ACTION(std::vector<float>, elePAtVtx, ## __VA_ARGS__)                \
   ACTION(std::vector<float>, elePAtSC, ## __VA_ARGS__)                 \
   ACTION(std::vector<float>, elePAtCluster, ## __VA_ARGS__)            \
   ACTION(std::vector<float>, elePAtSeed, ## __VA_ARGS__)               \
   ACTION(std::vector<float>, eleBrem, ## __VA_ARGS__)                  \
   ACTION(std::vector<float>, eledEtaAtVtx, ## __VA_ARGS__)             \
   ACTION(std::vector<float>, eledPhiAtVtx, ## __VA_ARGS__)             \
   ACTION(std::vector<float>, eleSigmaIEtaIEta, ## __VA_ARGS__)         \
   ACTION(std::vector<float>, eleSigmaIEtaIEta_2012, ## __VA_ARGS__)    \
   ACTION(std::vector<float>, eleSigmaIPhiIPhi, ## __VA_ARGS__)         \
   ACTION(std::vector<int>, eleMissHits, ## __VA_ARGS__)                \
   ACTION(std::vector<float>, eleESEffSigmaRR, ## __VA_ARGS__)          \
   ACTION(std::vector<float>, elePFChIso, ## __VA_ARGS__)               \
   ACTION(std::vector<float>, elePFPhoIso, ## __VA_ARGS__)              \
   ACTION(std::vector<float>, elePFNeuIso, ## __VA_ARGS__)              \
   ACTION(std::vector<float>, elePFPUIso, ## __VA_ARGS__)               \
   ACTION(std::vector<float>, elePFChIso03, ## __VA_ARGS__)             \
   ACTION(std::vector<float>, elePFPhoIso03, ## __VA_ARGS__)            \
   ACTION(std::vector<float>, elePFNeuIso03, ## __VA_ARGS__)            \
   ACTION(std::vector<float>, elePFChIso04, ## __VA_ARGS__)             \
   ACTION(std::vector<float>, elePFPhoIso04, ## __VA_ARGS__)            \
   ACTION(std::vector<float>, elePFNeuIso04, ## __VA_ARGS__)            \
   ACTION(std::vector<float>, eleR9, ## __VA_ARGS__)                    \
   ACTION(std::vector<float>, eleE3x3, ## __VA_ARGS__)                  \
   ACTION(std::vector<float>, eleE5x5, ## __VA_ARGS__)                  \
   ACTION(std::vector<float>, eleR9Full5x5, ## __VA_ARGS__)             \
   ACTION(std::vector<float>, eleE3x3Full5x5, ## __VA_ARGS__)           \
   ACTION(std::vector<float>, eleE5x5Full5x5, ## __VA_ARGS__)           \

#define B_VEC_L(ACTION, ...)                                            \
   ACTION(std::vector<float>, egEt, ## __VA_ARGS__)                     \
   ACTION(std::vector<float>, egEta, ## __VA_ARGS__)                    \
   ACTION(std::vector<float>, egPhi, ## __VA_ARGS__)                    \
   ACTION(std::vector<short>, egIEt, ## __VA_ARGS__)                    \
   ACTION(std::vector<short>, egIEta, ## __VA_ARGS__)                   \
   ACTION(std::vector<short>, egIPhi, ## __VA_ARGS__)                   \
   ACTION(std::vector<short>, egIso, ## __VA_ARGS__)                    \
   ACTION(std::vector<short>, egBx, ## __VA_ARGS__)                     \
   ACTION(std::vector<short>, egTowerIPhi, ## __VA_ARGS__)              \
   ACTION(std::vector<short>, egTowerIEta, ## __VA_ARGS__)              \
   ACTION(std::vector<short>, egRawEt, ## __VA_ARGS__)                  \
   ACTION(std::vector<short>, egIsoEt, ## __VA_ARGS__)                  \
   ACTION(std::vector<short>, egFootprintEt, ## __VA_ARGS__)            \
   ACTION(std::vector<short>, egNTT, ## __VA_ARGS__)                    \
   ACTION(std::vector<short>, egShape, ## __VA_ARGS__)                  \
   ACTION(std::vector<short>, egTowerHoE, ## __VA_ARGS__)               \
   ACTION(std::vector<short>, egHwQual, ## __VA_ARGS__)                 \

#define B_VEC_M(ACTION, ...)                                            \
   ACTION(std::vector<int>, mcPID, ## __VA_ARGS__)                      \
   ACTION(std::vector<int>, mcStatus, ## __VA_ARGS__)                   \
   ACTION(std::vector<float>, mcVtx_x, ## __VA_ARGS__)                  \
   ACTION(std::vector<float>, mcVtx_y, ## __VA_ARGS__)                  \
   ACTION(std::vector<float>, mcVtx_z, ## __VA_ARGS__)                  \
   ACTION(std::vector<float>, mcPt, ## __VA_ARGS__)                     \
   ACTION(std::vector<float>, mcEta, ## __VA_ARGS__)                    \
   ACTION(std::vector<float>, mcPhi, ## __VA_ARGS__)                    \
   ACTION(std::vector<float>, mcE, ## __VA_ARGS__)                      \
   ACTION(std::vector<float>, mcEt, ## __VA_ARGS__)                     \
   ACTION(std::vector<float>, mcMass, ## __VA_ARGS__)                   \
   ACTION(std::vector<int>, mcParentage, ## __VA_ARGS__)                \
   ACTION(std::vector<int>, mcMomPID, ## __VA_ARGS__)                   \
   ACTION(std::vector<float>, mcMomPt, ## __VA_ARGS__)                  \
   ACTION(std::vector<float>, mcMomEta, ## __VA_ARGS__)                 \
   ACTION(std::vector<float>, mcMomPhi, ## __VA_ARGS__)                 \
   ACTION(std::vector<float>, mcMomMass, ## __VA_ARGS__)                \
   ACTION(std::vector<int>, mcGMomPID, ## __VA_ARGS__)                  \
   ACTION(std::vector<int>, mcIndex, ## __VA_ARGS__)                    \
   ACTION(std::vector<float>, mcCalIsoDR03, ## __VA_ARGS__)             \
   ACTION(std::vector<float>, mcCalIsoDR04, ## __VA_ARGS__)             \
   ACTION(std::vector<float>, mcTrkIsoDR03, ## __VA_ARGS__)             \
   ACTION(std::vector<float>, mcTrkIsoDR04, ## __VA_ARGS__)             \

#define B_VEC_N(ACTION, ...)                                            \
   ACTION(std::vector<float>, eleTrkPtRelErr, ## __VA_ARGS__)           \
   ACTION(std::vector<int>, eleGenMatchIndex, ## __VA_ARGS__)           \
   ACTION(std::vector<int>, mcRecoMatchIndex, ## __VA_ARGS__)           \
   ACTION(std::vector<int>, hlt, ## __VA_ARGS__)                        \
   ACTION(std::vector<double>, pt_e20, ## __VA_ARGS__)                  \
   ACTION(std::vector<double>, eta_e20, ## __VA_ARGS__)                 \
   ACTION(std::vector<double>, phi_e20, ## __VA_ARGS__)                 \
   ACTION(std::vector<int>, n_e20, ## __VA_ARGS__)                      \

class electrontree {
   public:
      electrontree(bool do_mc_branches, bool do_l1_branches) {
         this->do_mc_branches = do_mc_branches;
         this->do_l1_branches = do_l1_branches;

         B_VAR_D(INVALID)
         B_VEC_D(NEWVEC)
         if (do_l1_branches) {
            B_VAR_L(INVALID)
            B_VEC_L(NEWVEC) }
         if (do_mc_branches) {
            B_VAR_M(INVALID)
            B_VEC_M(NEWVEC) }
         B_VAR_N(INVALID)
         B_VEC_N(NEWVEC)
      };

      electrontree(TTree* t, bool do_mc_branches, bool do_l1_branches)
            : electrontree(do_mc_branches, do_l1_branches) {
         branch(t);
      };

      electrontree(bool do_mc_branches, bool do_l1_branches, TTree* t)
            : electrontree(do_mc_branches, do_l1_branches) {
         read(t);
      };

      ~electrontree() {};

      void branch(TTree* t) {
         B_VAR_D(BRNREF, t)
         B_VEC_D(BRNVAR, t)
         if (do_l1_branches) {
            B_VAR_L(BRNREF, t)
            B_VEC_L(BRNVAR, t) }
         if (do_mc_branches) {
            B_VAR_M(BRNREF, t)
            B_VEC_M(BRNVAR, t) }
         B_VAR_N(BRNREF, t)
         B_VEC_N(BRNVAR, t)
      };

      void read(TTree* t) {
         B_VAR_D(RREF, t)
         B_VEC_D(RVAR, t)
         if (do_l1_branches) {
            B_VAR_L(RREF, t)
            B_VEC_L(RVAR, t) }
         if (do_mc_branches) {
            B_VAR_M(RREF, t)
            B_VEC_M(RVAR, t) }
         B_VAR_N(RREF, t)
         B_VEC_N(RVAR, t)
      };

      void clear() {
         B_VEC_D(CLEAR)
         if (do_l1_branches) {
            B_VEC_L(CLEAR) }
         if (do_mc_branches) {
            B_VEC_M(CLEAR) }
         B_VEC_N(CLEAR)
      };

      void copy(eventtree* evtt) {
         B_VAR_D(VARCOPY, evtt)
         B_VEC_D(VECCOPY, evtt)
         if (do_mc_branches) {
            B_VAR_M(VARCOPY, evtt)
            B_VEC_M(VECCOPY, evtt) }
      };

      void copy(l1tree* l1t) {
         if (do_l1_branches) {
            B_VAR_L(VARCOPY, l1t)
            B_VEC_L(VECCOPY, l1t)
         }
      };

      B_VAR_D(DECLARE)
      B_VAR_L(DECLARE)
      B_VAR_M(DECLARE)
      B_VAR_N(DECLARE)
      B_VEC_D(DECLPTR)
      B_VEC_L(DECLPTR)
      B_VEC_M(DECLPTR)
      B_VEC_N(DECLPTR)

   private:
      bool do_mc_branches;
      bool do_l1_branches;
};

#endif  /* ELECTRONTREE_H */
