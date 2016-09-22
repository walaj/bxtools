#include "bxstats.h"

#include "bxcommon.h"

#include <getopt.h>
#include <iostream>
#include <sstream>

#include "SeqLib/BamReader.h"

namespace opt {

  static std::string bam; // the bam to analyze
  static bool verbose = false; 

}

static const char* shortopts = "hv";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "bam",                     required_argument, NULL, 'b' },
  { NULL, 0, NULL, 0 }
};

static const char *STAT_USAGE_MESSAGE =
"Usage: bxtools stat <BAM/SAM/CRAM> > stats.tsv\n"
"Description: Gather BX-level statistics\n"
"\n"
"  General options\n"
"  -v, --verbose                        Set verbose output\n"
"\n";

void runStat(int argc, char** argv) {
  
  parseStatOptions(argc, argv);

  // open the BAM
  SeqLib::BamReader reader;
  if (!reader.Open(opt::bam)) {
    std::cerr << "Failed to open bam: " << opt::bam << std::endl;
    exit(EXIT_FAILURE);
  }

  std::unordered_map<std::string, BXStat> bxstats;

  // loop and collect
  SeqLib::BamRecord r;
  size_t count = 0;
  while (reader.GetNextRecord(r)) {
    
    ++count;

    // sanity check
    if (count == 100000 && !bxstats.size())
      std::cerr << "****1e5 reads in and haven't hit BX tag yet****" << std::endl;
    
    std::string bx = r.GetZTag("BX");
    if (bx.empty()) {
      //std::cerr << "BX tag empty for read: " << r << std::endl;
      continue;
    }
      
    if (count % 1000000 == 0)
      std::cerr << "...at read " << SeqLib::AddCommas(count) << " at pos " << r.Brief() << std::endl;
    
    ++bxstats[bx].count;
    bxstats[bx].bx = bx;
    if (r.PairMappedFlag() && !r.Interchromosomal())
      bxstats[bx].isize.push_back(std::abs(r.InsertSize()));
    if (r.MappedFlag())
      bxstats[bx].mapq.push_back(std::abs(r.MapQuality()));
    
  }

  for (const auto& b : bxstats)
    std::cout << b.second << std::endl;

}

void parseStatOptions(int argc, char** argv) {

  bool die = false;
  bool help = false;

  if (argc < 2) 
    die = true;

  opt::bam = std::string(argv[1]);

  std::stringstream ss;

  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v': opt::verbose = true; break;
    case 'h': help = true; break;
    }
  }

  if (die || help) 
    {
      std::cerr << "\n" << STAT_USAGE_MESSAGE;
      if (die)
	exit(EXIT_FAILURE);
      else 
	exit(EXIT_SUCCESS);	
    }

}

static double CalcMHWScore(std::vector<int> scores) {
  double median;
  size_t size = scores.size();

  std::sort(scores.begin(), scores.end());

  if (size  % 2 == 0)
    {
      median = (scores[size / 2 - 1] + scores[size / 2]) / 2;
    }
  else 
    {
      median = scores[size / 2];
    }

  return median;
}
