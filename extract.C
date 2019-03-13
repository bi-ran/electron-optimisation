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

TChain* chain_from_files(std::vector<std::string>& files, std::string name,
                         bool flag) {
   if (!flag) return nullptr;

   TChain* c = new TChain(name.data());
   for (auto const& f : files)
      c->Add(f.data());

   c->SetBranchStatus("*", 0);

   return c;
}

int extract(const char* config, const char* output) {
   configurer* conf = new configurer(config);

   auto files = conf->get<std::vector<std::string>>("files");
   auto paths = conf->get<std::vector<std::string>>("paths");
   auto mc_branches = conf->get<bool>("mc_branches");
   auto l1_branches = conf->get<bool>("l1_branches");
   auto hlt_branches = conf->get<bool>("hlt_branches");

   auto max_entries = conf->get<int64_t>("max_entries");

   TChain* ceg = chain_from_files(files, "ggHiNtuplizerGED/EventTree", true);
   TChain* cevt = chain_from_files(files, "hiEvtAnalyzer/HiTree", true);
   TChain* cl1 = chain_from_files(files, "l1object/L1UpgradeFlatTree", l1_branches);
   TChain* chlt = chain_from_files(files, "hltanalysis/HltTree", hlt_branches);
   TChain* cobj = chain_from_files(files, "hltobject/HLT_HIEle20Gsf_v", hlt_branches);

   int hiBin; RREF(int, hiBin, cevt);
   float hiHF; RREF(float, hiHF, cevt);

   std::vector<int> hlt(paths.size());
   if (hlt_branches) {
      for (std::size_t i=0; i<paths.size(); ++i) {
         chlt->SetBranchStatus(paths[i].data(), 1);
         chlt->SetBranchAddress(paths[i].data(), &hlt[i]);
      }
   }

   std::vector<double>* obj_pt = 0;
   std::vector<double>* obj_eta = 0;
   std::vector<double>* obj_phi = 0;

   if (hlt_branches) {
      cobj->SetBranchStatus("pt", 1);
      cobj->SetBranchAddress("pt", &obj_pt);
      cobj->SetBranchAddress("eta", &obj_eta);
      cobj->SetBranchAddress("phi", &obj_phi);
   }

   eventtree* evtt = new eventtree(ceg, mc_branches);
   l1tree* l1t = (l1_branches) ? new l1tree(cl1, l1_branches) : 0;

   TFile* fout = new TFile(output, "recreate");
   TTree* tout = new TTree("electrons", "electrons");

   electrontree* elet = new electrontree(tout, mc_branches, l1_branches);

   const float mindr2 = 0.15 * 0.15;

   int64_t nentries = ceg->GetEntries();
   if (max_entries > 0) { nentries = std::min(nentries, max_entries); }
   printf("entries: %lu\n", nentries);
   for (int64_t i=0; i<nentries; ++i) {
      elet->clear();

      ceg->GetEntry(i);
      cevt->GetEntry(i);

      if (hlt_branches) {
         chlt->GetEntry(i);
         cobj->GetEntry(i);
      }

      if (l1_branches)
         cl1->GetEntry(i);

      if (i % 10000 == 0) { printf("entry: %lu\n", i); }

      if (hlt_branches)
         for (auto const& b : hlt)
            elet->hlt->push_back(b);

      if (mc_branches) {
         elet->mcRecoMatchIndex->assign(evtt->nMC, -1);
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

            elet->eleGenMatchIndex->push_back(match);
            if (match != -1) { (*elet->mcRecoMatchIndex)[match] = j; }
         }
      }

      if (l1_branches)
         elet->copy(l1t);
      elet->copy(evtt);
      elet->hiBin = hiBin;
      elet->hiHF = hiHF;
      elet->ncoll = mc_branches ? ncoll(hiBin) : 1;

      for (int j=0; j<elet->nEle; ++j) {
         elet->eleTrkPtRelErr->push_back(
            (*elet->eleTrkPtErr)[j] / (*elet->eleTrkPt)[j]);
      }

      if (hlt_branches) {
         for (uint32_t j = 0; j < obj_pt->size(); ++j) {
            auto const& o = (*obj_pt)[i];
            if (std::find(elet->pt_obj->begin(), elet->pt_obj->end(), o)
                  == elet->pt_obj->end()) {
               elet->n_obj->push_back(std::count(
                  obj_pt->begin(), obj_pt->end(), o));
               elet->pt_obj->push_back((*obj_pt)[j]);
               elet->eta_obj->push_back((*obj_eta)[j]);
               elet->phi_obj->push_back((*obj_phi)[j]);
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
