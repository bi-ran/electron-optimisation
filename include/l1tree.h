#ifndef L1TREE_H
#define L1TREE_H

#include "TTree.h"

#include <vector>

#include "common.h"

#define B_AVL_L(ACTION, ...)                                            \
   ACTION(short, nEGs, ## __VA_ARGS__)                                  \

#define B_ARL_L(ACTION, ...)                                            \
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

class l1tree {
   public:
      l1tree(bool do_l1_branches) {
         this->do_l1_branches = do_l1_branches;

         B_AVL_L(ZERO)
         B_ARL_L(ZERO)
      };

      l1tree(TTree* t, bool do_l1_branches)
            : l1tree(do_l1_branches) {
         read(t);
      };

      ~l1tree() { };

      void read(TTree* t) {
         if (do_l1_branches) {
            B_AVL_L(RREF, t)
            B_ARL_L(RVAR, t) }
      };

      B_AVL_L(DECLARE)
      B_ARL_L(DECLPTR)

   private:
      bool do_l1_branches;
};

#endif  /* L1TREE_H */
