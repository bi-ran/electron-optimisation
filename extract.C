#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TLorentzVector.h"

#include "include/constants.h"
#include "include/electrontree.h"
#include "include/eventtree.h"
#include "include/l1tree.h"

#include "git/config/include/configurer.h"

#define PI 3.141593f

inline float dphi_2s1f1b(float phi1, float phi2) {
   float dphi = fabs(phi1 - phi2);
   if (dphi > PI) { dphi = 2 * PI - dphi; }
   return dphi;
}

int extract(const char* config, const char* output) {
   configurer* conf = new configurer(config);

   auto files = conf->get<std::vector<std::string>>("files");
   auto paths = conf->get<std::vector<std::string>>("paths");
   auto do_mc_branches = conf->get<bool>("do_mc_branches");
   auto do_l1_branches = conf->get<bool>("do_l1_branches");
   auto do_hlt_branches = conf->get<bool>("do_hlt_branches");

   auto maxentries = conf->get<uint64_t>("maxentries");

   TChain* ceg = new TChain("ggHiNtuplizerGED/EventTree");
   TChain* cevt = new TChain("hiEvtAnalyzer/HiTree");
   TChain* cl1 = (do_l1_branches) ? new TChain("l1object/L1UpgradeFlatTree") : 0;
   TChain* chlt = (do_hlt_branches) ? new TChain("hltanalysis/HltTree") : 0;
   TChain* ce20 = (do_hlt_branches) ? new TChain("hltobject/HLT_HIEle20Gsf_v") : 0;

   for (const auto& file : files) {
      ceg->Add(file.data());
      cevt->Add(file.data());

      if (do_l1_branches)
         cl1->Add(file.data());

      if (do_hlt_branches) {
         chlt->Add(file.data());
         ce20->Add(file.data());
      }
   }

   ceg->SetBranchStatus("*", 0);
   cevt->SetBranchStatus("*", 0);

   if (do_l1_branches)
      cl1->SetBranchStatus("*", 0);

   if (do_hlt_branches) {
      chlt->SetBranchStatus("*", 0);
      ce20->SetBranchStatus("*", 0);
   }

   int hiBin;
   cevt->SetBranchStatus("hiBin", 1);
   cevt->SetBranchAddress("hiBin", &hiBin);

   float hiHF;
   cevt->SetBranchStatus("hiHF", 1);
   cevt->SetBranchAddress("hiHF", &hiHF);

   std::vector<int> hlt(paths.size());
   if (do_hlt_branches) {
      for (std::size_t i=0; i<paths.size(); ++i) {
         chlt->SetBranchStatus(paths[i].data(), 1);
         chlt->SetBranchAddress(paths[i].data(), &hlt[i]);
      }
   }

   std::vector<double>* e20_pt = 0;
   std::vector<double>* e20_eta = 0;
   std::vector<double>* e20_phi = 0;

   if (do_hlt_branches) {
      ce20->SetBranchStatus("pt", 1);
      ce20->SetBranchAddress("pt", &e20_pt);
      ce20->SetBranchAddress("eta", &e20_eta);
      ce20->SetBranchAddress("phi", &e20_phi);
   }

   eventtree* evtt = new eventtree(ceg, do_mc_branches);
   l1tree* l1t = (do_l1_branches) ? new l1tree(cl1, do_l1_branches) : 0;

   TFile* fout = new TFile(output, "recreate");
   TTree* tout = new TTree("electrons", "electrons");

   electrontree* elet = new electrontree(tout, do_l1_branches, do_mc_branches);

   const float mindr2 = 0.15 * 0.15;

   uint64_t nentries = ceg->GetEntries();
   if (maxentries > 0) { nentries = std::min(nentries, maxentries); }
   printf("entries: %lu\n", nentries);
   for (uint64_t i=0; i<nentries; ++i) {
      elet->clear();

      ceg->GetEntry(i);
      cevt->GetEntry(i);

      if (do_hlt_branches) {
         chlt->GetEntry(i);
         ce20->GetEntry(i);
      }

      if (do_l1_branches)
         cl1->GetEntry(i);

      if (i % 10000 == 0) { printf("entry: %lu\n", i); }

      for (auto const& b : hlt) { elet->hlt.push_back(b); }

      if (do_mc_branches) {
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
      }

      if (do_l1_branches)
         elet->copy(l1t);
      elet->copy(evtt);
      elet->hiBin = hiBin;
      elet->hiHF = hiHF;
      elet->ncoll = do_mc_branches ? ncoll(hiBin) : 1;

      if (do_hlt_branches) {
         for (uint32_t j = 0; j < e20_pt->size(); ++j) {
            auto const& o = (*e20_pt)[i];
            if (std::find(elet->pt_e20.begin(), elet->pt_e20.end(), o)
                  == elet->pt_e20.end()) {
               elet->n_e20.push_back(std::count(
                  e20_pt->begin(), e20_pt->end(), o));
               elet->pt_e20.push_back((*e20_pt)[j]);
               elet->eta_e20.push_back((*e20_eta)[j]);
               elet->phi_e20.push_back((*e20_phi)[j]);
            }
         }
      }

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
