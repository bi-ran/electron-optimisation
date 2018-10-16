#include <fstream>
#include <string>

#include "git/config/configurer.h"

int generate(const char* config, const char* output) {
   configurer* conf = new configurer(config);

   auto title = conf->get<std::string>("title");
   auto overview = conf->get<std::vector<std::string>>("overview");
   auto figures = conf->get<std::vector<std::string>>("figures");
   auto labels = conf->get<std::vector<std::string>>("labels");
   auto captions = conf->get<std::vector<std::string>>("captions");
   auto perpage = conf->get<uint32_t>("perpage");

   auto nfigs = figures.size();

   if (!perpage || !nfigs) {
      printf("invalid values for perpage/figures!\n");
      return 1;
   }

   std::string outstr = "tex/electron-reco-comparison-";
   outstr += output; outstr += ".tex";

   std::fstream fs(outstr, std::fstream::out);
   fs << "\\documentclass[pdf]{beamer}\n"
         "\\mode<presentation>{\\usetheme{CambridgeUS}}\n"
         "\\setbeamertemplate{navigation symbols}{}\n"
         "\\setbeamertemplate{caption}{\\raggedright\\insertcaption\\par}\n"
         "\n"
         "\\title{Electron performance}\n"
         "\\subtitle{" << title << "}\n"
         "\n"
         "\\begin{document}\n"
         "\n"
         "\\begin{frame}\n"
         "\\titlepage\n"
         "\\end{frame}\n"
         "\n"
         "\\begin{frame}{General}\n"
         "\\begin{itemize}\n"
         "\\item Sample: $Z \\rightarrow e^{+}e^{-}$\n"
         "\\item Electron selections:\n"
         "    \\begin{itemize}\n"
         "    \\item electron $p_{T} > 10$\\,GeV\n"
         "    \\item 2015 veto ID cuts\n"
         "    \\end{itemize}\n";
   for (const auto& l : overview)
      fs << "\\item " << l << "\n";
   fs << "\\end{itemize}\n"
         "\\end{frame}\n"
         "\n";

   uint32_t f = 0;
   for (uint32_t i=0; i<(nfigs+perpage-1)/perpage; ++i) {
      fs << "\\begin{frame}{";
      for (uint32_t j=0, fl=f; j<perpage && fl<nfigs; ++j, ++fl) {
         fs << labels[fl]; if (j<perpage-1 && fl<nfigs-1) fs << ", "; }
      fs << "}\n"
            "\\begin{figure}[htb]\n"
            "\\centering\n";
      for (uint32_t j=0; j<perpage && f<nfigs; ++j, ++f)
         fs << "\\begin{minipage}{" << 0.99/perpage-0.01 << "\\textwidth}\n"
               "\\centering\n"
               "\\includegraphics[width=0.99\\textwidth]{" << figures[f] << "}\n"
               "\\caption{" << captions[f] << "}\n"
               "\\end{minipage}\n";
      fs << "\\end{figure}\n"
            "\\end{frame}\n"
            "\n";
   }

   fs << "\\end{document}\n";

   return 0;
}

int main(int argc, char* argv[]) {
   if (argc == 3)
      return generate(argv[1], argv[2]);

   printf("usage: %s [config] [label]\n", argv[0]);
   return 1;
}
