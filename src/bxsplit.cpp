#include "bxsplit.h"

#include <string>
#include <getopt.h>
#include <iostream>
#include <sstream>

namespace opt {

  static std::string bam; // the bam to split
  static std::string analysis_id; // unique prefix for output
}

static const char* shortopts = "hb:a:";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "bam",                     required_argument, NULL, 'b' },
  { "analysis-id",             required_argument, NULL, 'a' },
  { NULL, 0, NULL, 0 }
};

static const char *SPLIT_USAGE_MESSAGE =
"Usage: bxtools split -b <BAM/SAM/CRAM> -a <id>\n"
"Description: Split a BAM into multiple BAMs, one BAM per unique BX tag\n"
"\n"
"  General options\n"
"  -v, --verbose                        Select verbosity level (0-4). Default: 0 \n"
"  -h, --help                           Display this help and exit\n"
"\n";

void parseSplitOptions(int argc, char** argv) {

  bool die = false;

  if (argc <= 2) 
    die = true;

  bool help = false;
  std::stringstream ss;

  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'i': arg >> opt::bam; break;
    }

  if (die || help) 
    {
      std::cerr << "\n" << SPLIT_USAGE_MESSAGE;
      if (die)
	exit(EXIT_FAILURE);
      else 
	exit(EXIT_SUCCESS);	
    }

  }
}

void runSplit(int argc, char** argv) {
  
  parseSplitOptions(argc, argv);
  
  }
