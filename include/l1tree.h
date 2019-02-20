#ifndef L1TREE_H
#define L1TREE_H

#include "TTree.h"

#include <vector>

#define BRANCHES(ACTION)                              \
   BRANCHESDATA(ACTION)                               \

#define BRANCHESDATA(ACTION)                          \
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

#define ZERO(type, var) var = 0;
#define DECLARE(type, var) type var;
#define READ(type, var)                               \
   t->SetBranchStatus(#var, 1);                       \
   t->SetBranchAddress(#var, &var);                   \

class l1tree {
   public:
      l1tree() { BRANCHES(ZERO) };
      l1tree(TTree* t, bool isdata) : l1tree() { read(t, isdata); }
      ~l1tree() { };

      void read(TTree* t, bool isdata) {
         if (isdata) {
            BRANCHESDATA(READ) }
      };

      BRANCHES(DECLARE)
};

#undef BRANCHES
#undef BRANCHESMC
#undef BRANCHESDATA
#undef ZERO
#undef DECLARE
#undef READ

#endif  /* L1TREE_H */
