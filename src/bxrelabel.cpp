#include "bxrelabel.h"

#include <string>
#include <getopt.h>
#include <iostream>
#include <sstream>

#include "SeqLib/BamReader.h"
#include "SeqLib/BamWriter.h"

namespace opt {
  static std::string bam; // the bam to rename
  static bool verbose = false; 
}

static const char* shortopts = "hv";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "verbose",                 no_argument, NULL, 'v' },
  { NULL, 0, NULL, 0 }
};

static const char *RELABEL_USAGE_MESSAGE =
"Usage: bxtools relabel input.bam > relabeled.bam \n"
"Description: Move BX barcodes from BX tag to qname\n"
"\n"
"  General options\n"
"  -v, --verbose                        Select verbosity level (0-4). Default: 0 \n"
"  -h, --help                           Display this help and exit\n"
"\n";

void parseRelabelOptions(int argc, char** argv) {

  bool die = false;
  if (argc < 2) 
    die = true;
  else
    opt::bam = std::string(argv[1]);

  bool help = false;
  std::stringstream ss;

  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v': opt::verbose = true; break;
    }
  }

  if (die || help) {
    std::cerr << "\n" << RELABEL_USAGE_MESSAGE;
    die ? exit(EXIT_FAILURE) : exit(EXIT_SUCCESS);
  }
}

void runRelabel(int argc, char** argv) {

  parseRelabelOptions(argc, argv);
  
  // open the read BAM
  SeqLib::BamReader reader;
  if (!reader.Open(opt::bam)) {
    std::cerr << "Failed to open bam: " << opt::bam << std::endl;
    exit(EXIT_FAILURE);
  }

  // open the write BAM
  SeqLib::BamWriter w;
  if (!w.Open("-"))  {
    std::cerr << "Failed to open output stream" << std::endl;
    exit(EXIT_FAILURE);
  }
  w.SetHeader(reader.Header());
  w.WriteHeader();
  
  // loop and write
  SeqLib::BamRecord r;
  size_t count = 0;
  bool bxtaghit = false;
  while (reader.GetNextRecord(r)) {

    ++count;

    // sanity check
    if (count == 100000 && !bxtaghit)
      std::cerr << "****1e5 reads in and haven't hit BX tag yet****" << std::endl;

    std::string bx = r.GetZTag("BX");
    if (bx.empty()) {
      if (opt::verbose)
	std::cerr << "BX tag empty for read: " << r << std::endl;
      continue;
    } else {
      bxtaghit = true;
    }
      
    if (count % 1000000 == 0 && opt::verbose)
      std::cerr << "...at read " << SeqLib::AddCommas(count) << " at pos " << r.Brief() << std::endl;

    // set the read name with the BX tag, remove the old one
    r.SetQname(r.Qname() + "_" + bx);
    r.RemoveTag("BX");
    
    if (!w.WriteRecord(r)) {
      std::cerr << "failed to write read " << r << " to BAM for " << bx << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  
  w.Close();
}
