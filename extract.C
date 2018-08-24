#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TLorentzVector.h"

#include "include/eventtree.h"
#include "include/electrontree.h"

#include "git/config/configurer.h"

#define PI 3.141593f

inline float dphi_2s1f1b(float phi1, float phi2) {
   float dphi = fabs(phi1 - phi2);
   if (dphi > PI) { dphi = 2 * PI - dphi; }
   return dphi;
}

int extract(const char* config, const char* output) {
   configurer* conf = new configurer(config);
   auto files = conf->get<std::vector<std::string>>("files");

   TChain* ceg = new TChain("ggHiNtuplizerGED/EventTree");
   TChain* cevt = new TChain("hiEvtAnalyzer/HiTree");

   for (const auto& file : files) {
      ceg->Add(file.data());
      cevt->Add(file.data());
   }

   ceg->SetBranchStatus("*", 0);
   cevt->SetBranchStatus("*", 0);

   int hiBin;
   cevt->SetBranchStatus("hiBin", 1);
   cevt->SetBranchAddress("hiBin", &hiBin);
   float hiHF;
   cevt->SetBranchStatus("hiHF", 1);
   cevt->SetBranchAddress("hiHF", &hiHF);

   eventtree* evtt = new eventtree(ceg);

   TFile* fout = new TFile(output, "recreate");
   TTree* tout = new TTree("electrons", "electrons");

   electrontree* elet = new electrontree(tout);

   const float mindr2 = 0.15 * 0.15;

   uint64_t nentries = ceg->GetEntries();
   printf("total entries: %lu\n", nentries);
   for (uint64_t i=0; i<nentries; ++i) {
      elet->clear();

      ceg->GetEntry(i);
      cevt->GetEntry(i);
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

      float elePairZMass = -1;
      for (int j=0; j<evtt->nEle; ++j) {
         for (int k=j+1; k<evtt->nEle; ++k) {
            /* reconstruct Z peak */
            TLorentzVector e1; TLorentzVector e2;
            e1.SetPtEtaPhiM((*evtt->elePt)[j], (*evtt->eleEta)[j],
               (*evtt->elePhi)[j], 0.000511);
            e2.SetPtEtaPhiM((*evtt->elePt)[k], (*evtt->eleEta)[k],
               (*evtt->elePhi)[k], 0.000511);

            TLorentzVector zcand = e1 + e2;
            if (std::abs(zcand.M() - 91.1876) > elePairZMass)
               elePairZMass = zcand.M();
         }
      }

      elet->copy(evtt);
      elet->hiBin = hiBin;
      elet->hiHF = hiHF;
      elet->elePairZMass = elePairZMass;

      tout->Fill();
   }

   fout->Write("", TObject::kOverwrite);
   fout->Close();

   return 0;
}

int main(int argc, char* argv[]) {
   if (argc == 3)
      return extract(argv[1], argv[2]);

   printf("usage: %s [input] [output]\n", argv[0]);
   return 1;
}
