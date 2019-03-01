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
#include "TLine.h"

#include <algorithm>
#include <fstream>
#include <string>
#include <vector>

#include "include/cosmetics.h"

#include "git/config/include/configurer.h"

#define ASSERT(condition, message)     \
   if (!(condition)) {                 \
      printf(message "\n");            \
      return 1;                        \
   }

#define VECTOR_DEFAULT(var, n, val)          \
   if (var.empty()) { var.assign(n, val); }

int get_baseline(std::vector<uint32_t>& groups, uint32_t index);

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

   auto nbins = conf->get<std::vector<uint32_t>>("nbins");
   auto xbins = conf->get<std::vector<float>>("xbins");
   auto xrange = conf->get<std::vector<float>>("xrange");
   auto ybins = conf->get<std::vector<float>>("ybins");
   auto yrange = conf->get<std::vector<float>>("yrange");
   auto autorange = conf->get<bool>("autorange");

   auto csize = conf->get<std::vector<uint32_t>>("csize");
   auto logscale = conf->get<std::vector<bool>>("logscale");

   auto splitcanvas = conf->get<bool>("splitcanvas");
   auto drawratio = conf->get<bool>("drawratio");
   auto ratiotype = conf->get<int>("ratiotype");

   auto normalise = conf->get<int>("normalise");

   auto groups = conf->get<std::vector<uint32_t>>("groups");
   auto headers = conf->get<std::vector<std::string>>("headers");

   auto lines = conf->get<std::vector<float>>("lines");
   auto lncanvas = conf->get<std::vector<uint32_t>>("lncanvas");

   auto markers = conf->get<std::vector<int>>("markers");
   auto colours = conf->get<std::vector<int>>("colours");
   auto drawopts = conf->get<std::vector<std::string>>("drawopts");

   auto filename = conf->get<std::string>("filename");

   VECTOR_DEFAULT(xrange, 2, 0);
   VECTOR_DEFAULT(yrange, 2, 0);
   VECTOR_DEFAULT(logscale, 2, false);
   VECTOR_DEFAULT(drawopts, vars.size(), "p e same");

   ASSERT(!files.empty(), "no files provided")
   ASSERT(files.size() == trees.size() || trees.empty(),
      "#files != #trees")
   ASSERT(files.size() == vars.size() || vars.empty(),
      "#files != #vars")
   ASSERT(files.size() == selections.size() || selections.empty(),
      "#files != #selections")
   ASSERT(files.size() == tags.size(), "#files != #tags")
   ASSERT(files.size() == labels.size(), "#files != #labels")
   ASSERT(files.size() == legends.size(), "#files != #legends")
   ASSERT(xrange.size() == 2, "invalid x axis range")
   ASSERT(yrange.size() == 2, "invalid y axis range")
   ASSERT(csize.size() == 2, "invalid canvas size")
   ASSERT(groups.size() == headers.size(), "#groups != #headers")
   ASSERT(markers.size() == colours.size(), "#markers != #colours")
   ASSERT(colours.size() == files.size(), "#files != #colours")
   ASSERT(logscale.size() == 2, "(logx, logy)")
   ASSERT(lines.size() % 4 == 0, "line coordinates: n * (x0, y0, x1, y1)")

   if (!xbins.empty()) {
      ASSERT(xbins.size() == nbins[0] + 1, "#xbins != #nbins[0]+1") }

   if (drawratio && lines.empty()) {
      lncanvas.assign(1, splitcanvas);

      if (xbins.empty()) {
         lines.push_back(xrange[0]); lines.push_back(1);
         lines.push_back(xrange[1]); lines.push_back(1);
      } else {
         lines.push_back(xbins[0]); lines.push_back(1);
         lines.push_back(xbins[nbins[0]]); lines.push_back(1);
      }
   }

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

   TLine* ls[lines.size()/4];

   TH1* h[nfiles];

   for (std::size_t j = 0; j < nfiles; ++j) {
      f[j] = new TFile(files[j].data(), "read");
      t[j] = (TTree*)f[j]->Get(trees[j].data());

      switch (type) {
         case 0:
            if (xbins.empty())
               h1[j] = new TH1D(Form("hf%zu%s", j, tags[j].data()),
                     labels[j].data(), nbins[0], xrange[0], xrange[1]);
            else
               h1[j] = new TH1D(Form("hf%zu%s", j, tags[j].data()),
                     labels[j].data(), nbins[0], &xbins[0]);
            h[j] = h1[j];
            break;
         case 1:
            h2[j] = new TH2D(Form("hf%zu%s", j, tags[j].data()),
                  labels[j].data(), nbins[0], xrange[0], xrange[1],
                  nbins[1], yrange[0], yrange[1]);
            h[j] = h2[j];
            break;
         case 2:
            if (xbins.empty())
               hp[j] = new TProfile(Form("hf%zu%s", j, tags[j].data()),
                     labels[j].data(), nbins[0], xrange[0], xrange[1]);
            else
               hp[j] = new TProfile(Form("hf%zu%s", j, tags[j].data()),
                     labels[j].data(), nbins[0], &xbins[0]);
            h[j] = hp[j];
            break;
         case 3:
            h1[j] = (TH1D*)f[j]->Get(trees[j].data())->Clone(
               Form("hf%zu%s", j, tags[j].data()));
            break;
         default:
            break;
      }

      if (type < 3) {
         if (!common.empty()) selections[j] += (" && " + common);
         t[j]->Draw(Form("%s>>hf%zu%s", vars[j].data(), j, tags[j].data()),
               selections[j].data(), "goff");
      }

      switch (normalise) {
         case 0: break;
         case 1: h[j]->Scale(1. / h[j]->Integral()); break;
         case 2: h[j]->Scale(1. / h[j]->Integral(), "width"); break;
         case 3: h[j]->Scale(1. / t[j]->GetEntries()); break;
         case 4: h[j]->Scale(1. / t[j]->GetEntries(), "width"); break;
         case 5: h[j]->Scale(1. / t[j]->GetEntries(
            eventsel.data())); break;
         case 6: h[j]->Scale(1. / t[j]->GetEntries(
            eventsel.data()), "width"); break;
         default: break;
      }
   }

   hframe = (TH1D*)h[0]->Clone("hframe");

   float lmaxy = 0.84;
   float lminy = lmaxy - 0.045 * (nonempty + headers.size());
   TLegend* l1 = new TLegend(0.6, lminy, 0.96, lmaxy);
   lstyle(l1, 43, 12);

   for (std::size_t j = 0; j < nfiles; ++j) {
      hstyle(h[j], markers[j], colours[j], 0.7);

      if (nbins.size() == 1) {
         if (logscale[1]) { amin = std::min(amin, h[j]->GetMinimum(0)); }
         else { amin = std::min(amin, h[j]->GetMinimum()); }
         amax = std::max(amax, h[j]->GetMaximum());
      }

      unsigned k = std::abs(std::distance(groups.begin(), std::find(
         groups.begin(), groups.end(), j)));
      if (k < groups.size() && !headers[k].empty()) {
         TLegendEntry* e1 = l1->AddEntry((TObject*)0, headers[k].data(), "");
         e1->SetTextFont(63); e1->SetTextSize(12);
      }

      if (!legends[j].empty()) {
         l1->AddEntry(h[j], legends[j].data(), "pl"); }
   }

   if (autorange) {
      if (amax > amin * 1.08) amin = 0;
      hframe->SetAxisRange(amin * 0.8, amax * 1.44, "Y");
   } else if (!yrange.empty()) {
      hframe->SetAxisRange(yrange[0], yrange[1], "Y");
   }

   TCanvas* c1 = new TCanvas("c1", "", csize[0], csize[1]);
   if (drawratio && splitcanvas) {
      TPad* t1 = new TPad("p1", "", 0, 0.25, 1, 1);
      t1->SetTopMargin(0.11111); t1->SetBottomMargin(0);
      t1->Draw(); t1->SetNumber(1);
      TPad* t2 = new TPad("p2", "", 0, 0, 1, 0.25);
      t2->SetTopMargin(0); t2->SetBottomMargin(0.32);
      t2->Draw(); t2->SetNumber(2);

      hframe->GetXaxis()->SetLabelOffset(99);
      hframe->GetXaxis()->SetTitleOffset(99);
   }

   c1->cd(1);

   gPad->SetLogx(logscale[0]);

   if (!drawratio || splitcanvas) {
      gPad->SetLogy(logscale[1]);

      hstyle_title_label_size(hframe, 14, 12);
      hframe->Draw("axis");
   }

   if (drawratio) {
      c1->cd(splitcanvas + 1);

      hrframe = (TH1D*)h[0]->Clone("hrframe");
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

      hstyle_title_label_size(hrframe, 14, 12);
      hrframe->Draw("axis");
   }

   c1->cd(1);

   for (std::size_t j = 0; j < lines.size(); j=j+4) {
      if (j/4 < lncanvas.size()) { c1->cd(lncanvas[j/4] + 1); }

      ls[j/4] = new TLine(lines[j], lines[j+1], lines[j+2], lines[j+3]);
      ls[j/4]->SetLineStyle(7);
      ls[j/4]->Draw();
   }

   c1->cd(1);

   if (!drawratio || splitcanvas) {
      for (std::size_t j = 0; j < nfiles; ++j)
         h[j]->Draw(drawopts[j].data());
   }

   if (drawratio) {
      c1->cd(splitcanvas + 1);

      for (std::size_t j = 0; j < nfiles; ++j) {
         auto k = get_baseline(groups, j);
         if (k < 0 || k == (int)j) { continue; }

         switch (ratiotype) {
            case 0:
               hr1[j] = (TH1D*)h1[j]->Clone(
                  Form("hr1f%zu%s", j, tags[j].data()));
               hr1[j]->Divide(h1[k]);
               hr1[j]->Draw(drawopts[j].data());
               break;
            case 1:
               gr1[j] = new TGraphAsymmErrors(h1[j]->GetNbinsX() + 2);
               gr1[j]->SetName(Form("gr1f%zu%s", j, tags[j].data()));
               gr1[j]->Divide(h1[j], h1[k], "c1=0.683 b(1,1) mode");
               hstyle(gr1[j], markers[j], colours[j], 0.8);
               gr1[j]->Draw(drawopts[j].data());
            default:
               break;
         }
      }
   }

   c1->cd(1);

   l1->Draw();

   TLatex* t1 = new TLatex(); t1->SetTextFont(43); t1->SetTextSize(12);
   for (std::size_t l = 0; l < text.size(); ++l)
      t1->DrawLatexNDC(0.16, 0.825 - 0.03 * l, text[l].data());

   c1->SaveAs(Form("figs/%s-%s.pdf", filename.data(), output));
   c1->SaveAs(Form("figs/%s-%s.png", filename.data(), output));
   delete c1;

   TFile* fout = new TFile(Form("data/%s.root", output), "update");
   for (std::size_t j = 0; j < nfiles; ++j)
      h[j]->Write("", TObject::kOverwrite);
   fout->Close();

   return 0;
}

int get_baseline(std::vector<uint32_t>& groups, uint32_t index) {
   int max = -1;
   for (const auto& g : groups)
      if ((signed)(g - max) > 0 && g <= index)
         max = g;

   return max;
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
