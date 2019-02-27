#ifndef L1TREE_H
#define L1TREE_H

#include "TTree.h"

#include <vector>

#include "common.h"

#define BRANCHES_L1T(ACTION, ...)                                       \
   B_ALL_L(ACTION, ## __VA_ARGS__)                                      \

#define B_ALL_L(ACTION, ...)                                            \
   ACTION(short, nEGs, ## __VA_ARGS__)                                  \
   ACTION(std::vector<float>*, egEt, ## __VA_ARGS__)                    \
   ACTION(std::vector<float>*, egEta, ## __VA_ARGS__)                   \
   ACTION(std::vector<float>*, egPhi, ## __VA_ARGS__)                   \
   ACTION(std::vector<short>*, egIEt, ## __VA_ARGS__)                   \
   ACTION(std::vector<short>*, egIEta, ## __VA_ARGS__)                  \
   ACTION(std::vector<short>*, egIPhi, ## __VA_ARGS__)                  \
   ACTION(std::vector<short>*, egIso, ## __VA_ARGS__)                   \
   ACTION(std::vector<short>*, egBx, ## __VA_ARGS__)                    \
   ACTION(std::vector<short>*, egTowerIPhi, ## __VA_ARGS__)             \
   ACTION(std::vector<short>*, egTowerIEta, ## __VA_ARGS__)             \
   ACTION(std::vector<short>*, egRawEt, ## __VA_ARGS__)                 \
   ACTION(std::vector<short>*, egIsoEt, ## __VA_ARGS__)                 \
   ACTION(std::vector<short>*, egFootprintEt, ## __VA_ARGS__)           \
   ACTION(std::vector<short>*, egNTT, ## __VA_ARGS__)                   \
   ACTION(std::vector<short>*, egShape, ## __VA_ARGS__)                 \
   ACTION(std::vector<short>*, egTowerHoE, ## __VA_ARGS__)              \
   ACTION(std::vector<short>*, egHwQual, ## __VA_ARGS__)                \

class l1tree {
   public:
      l1tree() { do_l1_branches = 0; BRANCHES_L1T(ZERO) };
      l1tree(TTree* t, bool do_l1_branches)
         : l1tree() {
            this->do_l1_branches = do_l1_branches;
            read(t);
         }
      ~l1tree() { };

      void read(TTree* t) {
         if (do_l1_branches) {
            B_ALL_L(READ, t) }
      };

      bool do_l1_branches;

      BRANCHES_L1T(DECLARE)
};

#endif  /* L1TREE_H */
