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

constexpr int nbins = 200;
constexpr float Ncoll[nbins] = {1976.95, 1944.02, 1927.29, 1891.9, 1845.3, 1807.2, 1760.45, 1729.18, 1674.8, 1630.3, 1590.52, 1561.72, 1516.1, 1486.5, 1444.68, 1410.88, 1376.4, 1347.32, 1309.71, 1279.98, 1255.31, 1219.89, 1195.13, 1165.96, 1138.92, 1113.37, 1082.26, 1062.42, 1030.6, 1009.96, 980.229, 955.443, 936.501, 915.97, 892.063, 871.289, 847.364, 825.127, 806.584, 789.163, 765.42, 751.187, 733.001, 708.31, 690.972, 677.711, 660.682, 640.431, 623.839, 607.456, 593.307, 576.364, 560.967, 548.909, 530.475, 519.575, 505.105, 490.027, 478.133, 462.372, 451.115, 442.642, 425.76, 416.364, 405.154, 392.688, 380.565, 371.167, 360.28, 348.239, 340.587, 328.746, 320.268, 311.752, 300.742, 292.172, 281.361, 274.249, 267.025, 258.625, 249.931, 240.497, 235.423, 228.63, 219.854, 214.004, 205.425, 199.114, 193.618, 185.644, 180.923, 174.289, 169.641, 161.016, 157.398, 152.151, 147.425, 140.933, 135.924, 132.365, 127.017, 122.127, 117.817, 113.076, 109.055, 105.16, 101.323, 98.098, 95.0548, 90.729, 87.6495, 84.0899, 80.2237, 77.2201, 74.8848, 71.3554, 68.7745, 65.9911, 63.4136, 61.3859, 58.1903, 56.4155, 53.8486, 52.0196, 49.2921, 47.0735, 45.4345, 43.8434, 41.7181, 39.8988, 38.2262, 36.4435, 34.8984, 33.4664, 31.8056, 30.351, 29.2074, 27.6924, 26.7754, 25.4965, 24.2802, 22.9651, 22.0059, 21.0915, 19.9129, 19.1041, 18.1487, 17.3218, 16.5957, 15.5323, 14.8035, 14.2514, 13.3782, 12.8667, 12.2891, 11.61, 11.0026, 10.3747, 9.90294, 9.42648, 8.85324, 8.50121, 7.89834, 7.65197, 7.22768, 6.7755, 6.34855, 5.98336, 5.76555, 5.38056, 5.11024, 4.7748, 4.59117, 4.23247, 4.00814, 3.79607, 3.68702, 3.3767, 3.16309, 2.98282, 2.8095, 2.65875, 2.50561, 2.32516, 2.16357, 2.03235, 1.84061, 1.72628, 1.62305, 1.48916, 1.38784, 1.28366, 1.24693, 1.18552, 1.16085, 1.12596, 1.09298, 1.07402, 1.06105, 1.02954};

inline float ncoll(int hiBin) { return Ncoll[hiBin]; }

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
      elet->ncoll = ncoll(hiBin);
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
