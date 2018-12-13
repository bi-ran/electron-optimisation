#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TGraphAsymmErrors.h"
#include "TColor.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLatex.h"

#include <algorithm>
#include <fstream>
#include <string>
#include <vector>

#include "git/config/configurer.h"

#include "include/cosmetics.h"

#define ASSERT(condition, message)     \
   if (!(condition)) {                 \
      printf(message "\n");            \
      return 1;                        \
   }

#define VECTOR_DEFAULT(var, n, val)          \
   if (var.empty()) { var.assign(n, val); }

int get_baseline(std::vector<int> groups, uint32_t index);
void set_ratio_style(TH1D* h);

int harvest(const char* output, const char* config) {
   configurer* conf = new configurer(config);

   auto files = conf->get<std::vector<std::string>>("files");

   auto trees = conf->get<std::vector<std::string>>("trees");
   auto vars = conf->get<std::vector<std::string>>("vars");
   auto tags = conf->get<std::vector<std::string>>("tags");
   auto labels = conf->get<std::vector<std::string>>("labels");
   auto legends = conf->get<std::vector<std::string>>("legends");
   auto selections = conf->get<std::vector<std::string>>("selections");

   auto type = conf->get<int>("type");
   auto common = conf->get<std::string>("common");
   auto eventsel = conf->get<std::string>("eventsel");
   auto text = conf->get<std::vector<std::string>>("text");

   auto nbins = conf->get<std::vector<int>>("nbins");
   auto xbins = conf->get<std::vector<float>>("xbins");
   auto rmin = conf->get<std::vector<float>>("rmin");
   auto rmax = conf->get<std::vector<float>>("rmax");
   auto autorange = conf->get<bool>("autorange");

   auto csize = conf->get<std::vector<int>>("csize");
   auto logx = conf->get<bool>("logx");
   auto logy = conf->get<bool>("logy");

   auto splitcanvas = conf->get<bool>("splitcanvas");
   auto drawratio = conf->get<bool>("drawratio");
   auto ratiotype = conf->get<int>("ratiotype");

   auto normalise = conf->get<int>("normalise");

   auto groups = conf->get<std::vector<int>>("groups");
   auto headers = conf->get<std::vector<std::string>>("headers");

   auto markers = conf->get<std::vector<int>>("markers");
   auto colours = conf->get<std::vector<int>>("colours");

   auto filename = conf->get<std::string>("filename");

   ASSERT(!files.empty(), "no files provided")
   ASSERT(files.size() == trees.size(), "#files != #trees")
   ASSERT(files.size() == vars.size(), "#files != #vars")
   ASSERT(files.size() == tags.size(), "#files != #tags")
   ASSERT(files.size() == labels.size(), "#files != #labels")
   ASSERT(files.size() == legends.size(), "#files != #legends")
   ASSERT(files.size() == selections.size(), "#files != #selections")
   ASSERT(rmin.size() == 2, "invalid (min) ranges")
   ASSERT(rmax.size() == 2, "invalid (max) ranges")
   ASSERT(csize.size() == 2, "invalid canvas size")
   ASSERT(groups.size() == headers.size(), "#groups != #headers")
   ASSERT(markers.size() == colours.size(), "#markers != #colours")
   ASSERT(colours.size() == files.size(), "#files != #colours")

   std::size_t nfiles = files.size();

   std::for_each(legends.begin(), legends.end(), &rtrim);
   uint32_t nonempty = std::count_if(legends.begin(), legends.end(),
      [](std::string& l) { return !l.empty(); });

   double amin = 1;
   double amax = 0;

   TH1::SetDefaultSumw2();
   gStyle->SetOptStat(0);

   TFile* f[nfiles]; TTree* t[nfiles];

   TH1D* hframe; TH1D* hrframe;
   TH1D* h1[nfiles]; TH2D* h2[nfiles]; TProfile* hp[nfiles];
   TH1D* hr1[nfiles]; TGraphAsymmErrors* gr1[nfiles];

   TH1* h[nfiles];

   for (std::size_t j = 0; j < nfiles; ++j) {
      f[j] = new TFile(files[j].c_str(), "read");
      t[j] = (TTree*)f[j]->Get(trees[j].c_str());

      switch (type) {
         case 0:
            if (xbins.empty())
               h1[j] = new TH1D(Form("hf%zu%s", j, tags[j].c_str()),
                     labels[j].c_str(), nbins[0], rmin[0], rmax[0]);
            else
               h1[j] = new TH1D(Form("hf%zu%s", j, tags[j].c_str()),
                     labels[j].c_str(), nbins[0], &xbins[0]);
            h[j] = h1[j];
            break;
         case 1:
            h2[j] = new TH2D(Form("hf%zu%s", j, tags[j].c_str()),
                  labels[j].c_str(), nbins[0], rmin[0], rmax[0],
                  nbins[1], rmin[1], rmin[1]);
            h[j] = h2[j];
            break;
         case 2:
            if (xbins.empty())
               hp[j] = new TProfile(Form("hf%zu%s", j, tags[j].c_str()),
                     labels[j].c_str(), nbins[0], rmin[0], rmax[0]);
            else
               hp[j] = new TProfile(Form("hf%zu%s", j, tags[j].c_str()),
                     labels[j].c_str(), nbins[0], &xbins[0]);
            h[j] = hp[j];
            break;
         default:
            break;
      }

      if (!common.empty()) selections[j] += (" && " + common);
      t[j]->Draw(Form("%s>>hf%zu%s", vars[j].c_str(), j, tags[j].c_str()),
            selections[j].c_str(), "goff");

      switch (normalise) {
         case 0: break;
         case 1: h[j]->Scale(1. / h[j]->Integral()); break;
         case 2: h[j]->Scale(1. / h[j]->Integral(), "width"); break;
         case 3: h[j]->Scale(1. / t[j]->GetEntries()); break;
         case 4: h[j]->Scale(1. / t[j]->GetEntries(), "width"); break;
         case 5: h[j]->Scale(1. / t[j]->GetEntries(
            eventsel.c_str())); break;
         case 6: h[j]->Scale(1. / t[j]->GetEntries(
            eventsel.c_str()), "width"); break;
         default: break;
      }
   }

   hframe = (TH1D*)h[0]->Clone("hframe");

   float lmaxy = 0.835;
   float lminy = lmaxy - 0.04 * (nonempty + headers.size());
   TLegend* l1 = new TLegend(0.6, lminy, 0.96, lmaxy);
   lstyle(l1, 43, 12);

   for (std::size_t j = 0; j < nfiles; ++j) {
      hstyle(h[j], markers[j], colours[j], 0.8);

      if (nbins.size() == 1) {
         if (logy) { amin = std::min(amin, h[j]->GetMinimum(0)); }
         else { amin = std::min(amin, h[j]->GetMinimum()); }
         amax = std::max(amax, h[j]->GetMaximum());
      }

      unsigned k = std::abs(std::distance(groups.begin(), std::find(
         groups.begin(), groups.end(), j)));
      if (k < groups.size() && !headers[k].empty()) {
         TLegendEntry* e1 = l1->AddEntry((TObject*)0, headers[k].c_str(), "");
         e1->SetTextFont(63); e1->SetTextSize(13);
      }

      if (!legends[j].empty()) {
         l1->AddEntry(h[j], legends[j].c_str(), "pl"); }
   }

   if (autorange) {
      if (amax > amin * 1.08) amin = 0;
      hframe->SetAxisRange(amin * 0.8, amax * 1.44, "Y");
   } else if (rmin.size() > 1 && rmax.size() > 1) {
      hframe->SetAxisRange(rmin[1], rmax[1], "Y");
   }

   int cheight = splitcanvas ? csize[1] * 1.2 : csize[1];
   TCanvas* c1 = new TCanvas("c1", "", csize[0], cheight);
   if (drawratio && splitcanvas) {
      TPad* t1 = new TPad("p1", "", 0, 0.25, 1, 1);
      t1->SetTopMargin(0.11111); t1->SetBottomMargin(0);
      t1->Draw(); t1->SetNumber(1);
      TPad* t2 = new TPad("p2", "", 0, 0, 1, 0.25);
      t2->SetTopMargin(0); t2->SetBottomMargin(0.32);
      t2->Draw(); t2->SetNumber(2);
      c1->cd(1);

      hframe->GetXaxis()->SetLabelOffset(99);
      hframe->GetXaxis()->SetTitleOffset(99);
   }

   gPad->SetLogx(logx);
   gPad->SetLogy(logy);

   if (!drawratio || splitcanvas) {
      hframe->Draw("axis");
      for (std::size_t j = 0; j < nfiles; ++j)
         h[j]->Draw("p e same");
   }

   if (drawratio) {
      c1->cd(splitcanvas);

      hrframe = (TH1D*)h[0]->Clone("hrframe");
      set_ratio_style(hrframe);
      switch (ratiotype) {
         case 0:
            hrframe->SetAxisRange(0, 2, "Y");
            hrframe->SetNdivisions(205, "Y");
            hrframe->SetYTitle("ratio");
            break;
         case 1:
            hrframe->SetAxisRange(0, 1.2, "Y");
            hrframe->SetYTitle("efficiency");
            break;
      }
      hrframe->Draw("axis");

      for (std::size_t j = 0; j < nfiles; ++j) {
         auto k = get_baseline(groups, j);
         if (k < 0 || k == (int)j) { continue; }

         switch (ratiotype) {
            case 0:
               hr1[j] = (TH1D*)h1[j]->Clone(
                  Form("hr1f%zu%s", j, tags[j].c_str()));
               hr1[j]->Divide(h1[k]);
               hr1[j]->Draw("p e same");
               break;
            case 1:
               gr1[j] = new TGraphAsymmErrors(h1[j]->GetNbinsX() + 2);
               gr1[j]->SetName(Form("gr1f%zu%s", j, tags[j].c_str()));
               gr1[j]->Divide(h1[j], h1[k], "c1=0.683 b(1,1) mode");
               hstyle(gr1[j], markers[j], colours[j], 0.8);
               gr1[j]->Draw("p e same");
            default:
               break;
         }
      }
   }

   c1->cd(0);

   l1->Draw();

   TLatex* t1 = new TLatex(); t1->SetTextFont(43); t1->SetTextSize(13);
   for (std::size_t l = 0; l < text.size(); ++l)
      t1->DrawLatexNDC(0.16, 0.825 - 0.03 * l, text[l].c_str());

   c1->SaveAs(Form("figs/%s-%s.pdf", filename.c_str(), output));
   c1->SaveAs(Form("figs/%s-%s.png", filename.c_str(), output));
   delete c1;

   TFile* fout = new TFile(Form("data/%s.root", output), "update");
   for (std::size_t j = 0; j < nfiles; ++j)
      h[j]->Write("", TObject::kOverwrite);
   fout->Close();

   return 0;
}

int get_baseline(std::vector<int> groups, uint32_t index) {
   int max = -1;
   for (const auto& g : groups)
      if (g > max && g <= (int)index)
         max = g;

   return max;
}

void set_ratio_style(TH1D* h) {
   float npixelspad = gPad->GetWh() * gPad->GetAbsHNDC();

   h->GetXaxis()->SetLabelSize(13 / npixelspad);
   h->GetXaxis()->SetTitleSize(16 / npixelspad);
   h->GetYaxis()->SetLabelSize(13 / npixelspad);
   h->GetYaxis()->SetTitleSize(16 / npixelspad);
   h->GetYaxis()->CenterTitle();
   h->GetYaxis()->SetTitleOffset(0.5);
}

int main(int argc, char* argv[]) {
   if (argc > 2) {
      for (int f = 2; f < argc; ++f)
         harvest(argv[1], argv[f]);
   } else {
      printf("usage: %s [output] [configs ...]\n", argv[0]);
      return 1;
   }
}
