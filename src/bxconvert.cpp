#include "bxconvert.h"

#include <getopt.h>
#include <iostream>
#include <sstream>

#include "SeqLib/BamReader.h"
#include "SeqLib/BamWriter.h"
#include "SeqLib/GenomicRegionCollection.h"

#include "bxcommon.h"

namespace opt {

  static std::string bam; // the bam to analyze                                 
  static bool verbose = false;
  static int width = 1000;
  static int overlap = 0;
  static std::string bed; // optional bed file                                  
}

static const char* shortopts = "hvw:O:b:";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "bed",                     required_argument, NULL, 'b' },
  { "pad",                     required_argument, NULL, 'p' },
  { "width",                   required_argument, NULL, 'w' },
  { "overlap",                 required_argument, NULL, 'O' },
  { NULL, 0, NULL, 0 }
};


  void runConvert(int argc, char** argv) {

    parseConvertOptions(argc, argv);

    SeqLib::BamReader reader;
    BXOPEN(reader, opt::bam);
    SeqLib::BamHeader hdr = reader.Header();
    
    SeqLib::BamRecord r;
    size_t count = 0;
    std::unordered_map<std::string, size_t> bxtags;

    while (reader.GetNextRecord(r)){
      if (!bxtags.count(r.GetZTag("BX"))) {
        bxtags.insert(std::pair<std::string, size_t>(r.GetZTag("BX"), ++count));
      }
      std::string chr = hdr.IDtoName(r.ChrID());
      r.SetChrID(bxtags[r.GetZTag("BX")]);
      r.AddZTag("CR", chr);
        }

  }


void parseConvertOptions(int argc, char** argv) {

  bool die = false;
  bool help = false;

  if (argc < 2)
    die = true;

  opt::bam = std::string(argv[1]);

  std::stringstream ss;

  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)	\
    {
      std::istringstream arg(optarg != NULL ? optarg : "");
      switch (c) {
      case 'v': opt::verbose = true; break;
      case 'h': help = true; break;
      case 'w': arg >> opt::width; break;
      case 'O': arg >> opt::overlap; break;
      case 'b': arg >> opt::bed; break;
      }
    }

}
