#include "TFile.h"
#include "TTree.h"

#include "eventtree.h"
#include "electrontree.h"

#define PI 3.141593f

inline float dphi_2s1f1b(float phi1, float phi2) {
   float dphi = fabs(phi1 - phi2);
   if (dphi > PI) { dphi = 2 * PI - dphi; }
   return dphi;
}

int skim_electron(const char* input, const char* output) {
   TFile* fin = new TFile(input, "read");
   TTree* tin = (TTree*)fin->Get("ggHiNtuplizerGED/EventTree");
   tin->SetBranchStatus("*", 0);

   eventtree* evtt = new eventtree(tin);

   TFile* fout = new TFile(output, "recreate");
   TTree* tout = new TTree("electrons", "electrons");

   electrontree* elet = new electrontree(tout);

   const float mindr2 = 0.15 * 0.15;

   uint64_t nentries = tin->GetEntries();
   printf("total entries: %lu\n", nentries);
   for (uint64_t i=0; i<nentries; ++i) {
      elet->clear();

      tin->GetEntry(i);
      if (i % 10000 == 0) { printf("entry: %lu\n", i); }

      elet->mcRecoMatchIndex.assign(evtt->nMC, -1);

      for (int j=0; j<evtt->nEle; ++j) {
         int match = -1;
         float maxpt = -1;

         for (int k=0; k<evtt->nMC; ++k) {
            if (abs((*evtt->mcPID)[k]) != 11) { continue; }
            if ((*evtt->mcStatus)[k] != 1) { continue; }

            if ((*evtt->mcPt)[k] < maxpt) { continue; }

            float dphi = dphi_2s1f1b((*evtt->elePhi)[j], (*evtt->mcPhi)[k]);
            float deta = fabs((*evtt->eleEta)[j] - (*evtt->mcEta)[k]);
            float dr2 = dphi * dphi + deta * deta;

            if (dr2 < mindr2) {
               match = k;
               maxpt = (*evtt->mcPt)[k];
            }
         }

         elet->eleGenMatchIndex.push_back(match);
         if (match != -1) { elet->mcRecoMatchIndex[match] = j; }
      }

      elet->copy(evtt);

      tout->Fill();
   }

   fout->Write("", TObject::kOverwrite);
   fout->Close();

   return 0;
}

int main(int argc, char* argv[]) {
   if (argc == 3) {
      return skim_electron(argv[1], argv[2]);
   } else {
      printf("usage: ./skim_electron [input] [output]\n");
      return 1;
   }
}
