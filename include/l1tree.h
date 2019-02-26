#ifndef L1TREE_H
#define L1TREE_H

#include "TTree.h"

#include <vector>

#define L1BRANCHES(ACTION)                            \
   L1DATABRANCHES(ACTION)                             \

#define L1DATABRANCHES(ACTION)                        \
   ACTION(short, nEGs)                                \
   ACTION(std::vector<float>*, egEt)                  \
   ACTION(std::vector<float>*, egEta)                 \
   ACTION(std::vector<float>*, egPhi)                 \
   ACTION(std::vector<short>*, egIEt)                 \
   ACTION(std::vector<short>*, egIEta)                \
   ACTION(std::vector<short>*, egIPhi)                \
   ACTION(std::vector<short>*, egIso)                 \
   ACTION(std::vector<short>*, egBx)                  \
   ACTION(std::vector<short>*, egTowerIPhi)           \
   ACTION(std::vector<short>*, egTowerIEta)           \
   ACTION(std::vector<short>*, egRawEt)               \
   ACTION(std::vector<short>*, egIsoEt)               \
   ACTION(std::vector<short>*, egFootprintEt)         \
   ACTION(std::vector<short>*, egNTT)                 \
   ACTION(std::vector<short>*, egShape)               \
   ACTION(std::vector<short>*, egTowerHoE)            \
   ACTION(std::vector<short>*, egHwQual)              \

#define L1ZERO(type, var) var = 0;
#define L1DECLARE(type, var) type var;
#define L1READ(type, var)                             \
   t->SetBranchStatus(#var, 1);                       \
   t->SetBranchAddress(#var, &var);                   \

class l1tree {
   public:
      l1tree() { L1BRANCHES(L1ZERO) };
      l1tree(TTree* t, bool hasl1) : l1tree() { read(t, hasl1); }
      ~l1tree() { };

      void read(TTree* t, bool hasl1) {
         if (hasl1) {
            L1DATABRANCHES(L1READ) }
      };

      L1BRANCHES(L1DECLARE)
};

#endif  /* L1TREE_H */
