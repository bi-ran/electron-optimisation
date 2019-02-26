#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TLorentzVector.h"

#include "include/constants.h"
#include "include/eventtree.h"
#include "include/electrontree.h"
#include "include/l1tree.h"

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
   auto paths = conf->get<std::vector<std::string>>("paths");
   auto isdata = conf->get<bool>("isdata");

   auto maxentries = conf->get<uint64_t>("maxentries");

   TChain* ceg = new TChain("ggHiNtuplizerGED/EventTree");
   TChain* cevt = new TChain("hiEvtAnalyzer/HiTree");
   TChain* chlt = new TChain("hltanalysis/HltTree");
   TChain* ce20 = new TChain("hltobject/HLT_HIEle20Gsf_v");
   TChain* ce10e10m50 = new TChain("hltobject/HLT_HIDoubleEle10GsfMass50_v");
   TChain* cl1 = new TChain("l1object/L1UpgradeFlatTree");

   for (const auto& file : files) {
      ceg->Add(file.data());
      cevt->Add(file.data());
      chlt->Add(file.data());
      ce20->Add(file.data());
      ce10e10m50->Add(file.data());
      cl1->Add(file.data());
   }

   ceg->SetBranchStatus("*", 0);
   cevt->SetBranchStatus("*", 0);
   chlt->SetBranchStatus("*", 0);
   ce20->SetBranchStatus("*", 0);
   ce10e10m50->SetBranchStatus("*", 0);
   cl1->SetBranchStatus("*", 0);

   int hiBin;
   cevt->SetBranchStatus("hiBin", 1);
   cevt->SetBranchAddress("hiBin", &hiBin);

   float hiHF;
   cevt->SetBranchStatus("hiHF", 1);
   cevt->SetBranchAddress("hiHF", &hiHF);

   std::vector<int> hlt(paths.size());
   for (std::size_t i=0; i<paths.size(); ++i) {
      chlt->SetBranchStatus(paths[i].data(), 1);
      chlt->SetBranchAddress(paths[i].data(), &hlt[i]);
   }

   std::vector<double>* e20_pt = 0;
   std::vector<double>* e20_eta = 0;
   std::vector<double>* e20_phi = 0;
   ce20->SetBranchStatus("pt", 1);
   ce20->SetBranchAddress("pt", &e20_pt);
   ce20->SetBranchAddress("eta", &e20_eta);
   ce20->SetBranchAddress("phi", &e20_phi);

   std::vector<double>* e10e10m50_pt = 0;
   std::vector<double>* e10e10m50_eta = 0;
   std::vector<double>* e10e10m50_phi = 0;
   ce10e10m50->SetBranchStatus("pt", 1);
   ce10e10m50->SetBranchAddress("pt", &e10e10m50_pt);
   ce10e10m50->SetBranchAddress("eta", &e10e10m50_eta);
   ce10e10m50->SetBranchAddress("phi", &e10e10m50_phi);

   eventtree* evtt = new eventtree(ceg, isdata);
   l1tree* l1ot = new l1tree(cl1, isdata);

   TFile* fout = new TFile(output, "recreate");
   TTree* tout = new TTree("electrons", "electrons");

   electrontree* elet = new electrontree(tout, isdata);

   const float mindr2 = 0.15 * 0.15;

   uint64_t nentries = ceg->GetEntries();
   if (maxentries > 0) { nentries = std::min(nentries, maxentries); }
   printf("entries: %lu\n", nentries);
   for (uint64_t i=0; i<nentries; ++i) {
      elet->clear();

      ceg->GetEntry(i);
      cevt->GetEntry(i);
      chlt->GetEntry(i);
      ce20->GetEntry(i);
      ce10e10m50->GetEntry(i);
      cl1->GetEntry(i);

      if (i % 10000 == 0) { printf("entry: %lu\n", i); }

      for (auto const& b : hlt) { elet->hlt.push_back(b); }

      if (!elet->isData) {
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

      elet->copy(l1ot);
      elet->copy(evtt);
      elet->hiBin = hiBin;
      elet->hiHF = hiHF;
      elet->ncoll = elet->isData ? 1 : ncoll(hiBin);

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

      for (uint32_t j = 0; j < e10e10m50_pt->size(); ++j) {
         auto const& o = (*e10e10m50_pt)[i];
         if (std::find(elet->pt_e10e10m50.begin(), elet->pt_e10e10m50.end(), o)
               == elet->pt_e10e10m50.end()) {
            elet->n_e10e10m50.push_back(std::count(
               e10e10m50_pt->begin(), e10e10m50_pt->end(), o));
            elet->pt_e10e10m50.push_back((*e10e10m50_pt)[j]);
            elet->eta_e10e10m50.push_back((*e10e10m50_eta)[j]);
            elet->phi_e10e10m50.push_back((*e10e10m50_phi)[j]);
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
