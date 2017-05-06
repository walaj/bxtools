#include "bxgroup.h"

#include "bxcommon.h"
#include <string>
#include <getopt.h>
#include <iostream>
#include <sstream>

#include "SeqLib/BamReader.h"
#include "SeqLib/BamWriter.h"

struct BXGroup {

  int start;
  int stop;
  std::string tag;
  SeqLib::BamRecordVector buff;

  BXGroup() : start(0), stop(0), tag("MI") {}

  int width() const { return stop - start; }
};

namespace opt {

  static std::string bam;        // the bam to group
  static bool verbose = false; 
  static std::string tag = "BX"; // tag to group by
}

static const char* shortopts = "hvxb:a:m:t:";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "no-output",               no_argument, NULL, 'x' },
  { "verbose",                 no_argument, NULL, 'v' },
  { "min-reads",               required_argument, NULL, 'm' },
  { "tag",                     required_argument, NULL, 't' },
  { NULL, 0, NULL, 0 }
};

static const char *GROUP_USAGE_MESSAGE =
"Usage: bxtools group <BAM/SAM/CRAM> \n"
"Description: Group BX tags that are adjacent\n"
"\n"
"  General options\n"
"  -v, --verbose                        Select verbosity level (0-4). Default: 0 \n"
"  -h, --help                           Display this help and exit\n"
"  -a, --analysis-id                    ID to prefix output files with [foo]\n"
"  -x, --no-output                      Don't output BAMs (count only) [off]\n"
"  -m, --min-reads                      Minumum reads of given tag to see before writing [0]\n"
"  -t, --tag                            Split by a tag other than BX (e.g. MI)\n"
"\n";

void parseGroupOptions(int argc, char** argv) {

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
    case 'a': arg >> opt::analysis_id; break;
    case 'v': opt::verbose = true; break;
    case 't': arg >> opt::tag; break;
    }
  }

  if (die || help) {
    std::cerr << "\n" << GROUP_USAGE_MESSAGE;
    die ? exit(EXIT_FAILURE) : exit(EXIT_SUCCESS);
  }
}

void runGroup(int argc, char** argv) {
  
  parseGroupOptions(argc, argv);
  
  // opeen the BAM
  SeqLib::BamReader reader;
  if (!reader.Open(opt::bam)) {
    std::cerr << "Failed to open bam: " << opt::bam << std::endl;
    exit(EXIT_FAILURE);
  }
  
  // loop and write
  SeqLib::BamRecord r;
  size_t count = 0;
  bool hit = false;
  while (reader.GetNextRecord(r)) {

    ++count;

    // sanity check
    BXLOOPCHECK(r, hit, opt::tag)

    std::string bx;
    r.GetTag(opt::tag, bx);
    if (bx.empty()) {
      continue;
    } else {
      hit = true;
    }
    
  }  
}
