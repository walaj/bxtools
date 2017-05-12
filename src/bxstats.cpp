#include "bxstats.h"

#include "bxcommon.h"

#include <getopt.h>
#include <iostream>
#include <sstream>

#include "SeqLib/BamReader.h"

namespace opt {

  static std::string bam; // the bam to analyze
  static bool verbose = false; 
  static std::string tag = "BX"; // tag to split by
}

static const char* shortopts = "hvt:";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "tag",                     required_argument, NULL, 't' },
  { "bam",                     required_argument, NULL, 'b' },
  { NULL, 0, NULL, 0 }
};

static const char *STAT_USAGE_MESSAGE =
"Usage: bxtools stat <BAM/SAM/CRAM> > stats.tsv\n"
"Description: Gather BX-level statistics\n"
"\n"
"  General options\n"
"  -v, --verbose                        Set verbose output\n"
"  -t, --tag                            Collect stats by a tag other than BX (e.g. MI)\n"
"\n";

static void parseOptions(int argc, char** argv);

void runStat(int argc, char** argv) {
  
  parseOptions(argc, argv);

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
    std::string bx;
    bool tag_present = r.GetZTag(opt::tag, bx);
    BXLOOPCHECK(r, bxstats.size(), opt::tag)
    if (!tag_present)
      continue;

    ++bxstats[bx].count;
    bxstats[bx].bx = bx;
    if (r.PairMappedFlag() && !r.Interchromosomal())
      bxstats[bx].isize.push_back(std::abs(r.InsertSize()));
    if (r.MappedFlag())
      bxstats[bx].mapq.push_back(std::abs(r.MapQuality()));

    int as_int = -1;
    float as_float = -1;
    std::string as_string = "NA";
    if (r.GetIntTag("AS", as_int))
      bxstats[bx].as.push_back(as_int);
    else if (r.GetFloatTag("AS", as_float))
      bxstats[bx].as.push_back(as_float);
    else if (r.GetZTag("AS", as_string)) {
      try {
	bxstats[bx].as.push_back(std::stof(as_string));
      } catch (...) {
	std::cerr << "Could not convert AS:Z val of " << as_string << " to float" << std::endl;
      }
    }

  }

  for (const auto& b : bxstats)
    std::cout << b.second << std::endl;

}

static void parseOptions(int argc, char** argv) {

  bool die = false;
  bool help = false;

  if (argc < 2) 
    die = true;
  else
    opt::bam = std::string(argv[1]);

  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v': opt::verbose = true; break;
    case 'h': help = true; break;
    }
  }

  if (die || help) {
    std::cerr << "\n" << STAT_USAGE_MESSAGE;
    die ? exit(EXIT_FAILURE) : exit(EXIT_SUCCESS);	
  }

}

//http://stackoverflow.com/questions/2114797/compute-median-of-values-stored-in-vector-c
template <class T>
static double CalcMHWScore(std::vector<T> scores) {
  double median;
  size_t size = scores.size();

  std::sort(scores.begin(), scores.end());

  if (size  % 2 == 0)
      median = (scores[size / 2 - 1] + scores[size / 2]) / 2;
  else 
      median = scores[size / 2];

  return median;
}

std::ostream& operator<<(std::ostream& out, const BXStat& b) {
  double isize_med = -1;
  double mapq_med = -1;
  double as_med = -1;
  if (b.isize.size())
    isize_med = CalcMHWScore(b.isize);
  if (b.mapq.size())
    mapq_med = CalcMHWScore(b.mapq);
  if (b.as.size())
    as_med = CalcMHWScore(b.as);
  out << b.bx << "\t" << b.count << "\t" << isize_med << "\t" << mapq_med 
      << "\t" << as_med;
  return out;
}

