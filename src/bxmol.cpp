#include "bxmol.h"

#include <getopt.h>
#include <iostream>
#include <sstream>

#include "SeqLib/BamReader.h"
#include "SeqLib/GenomicRegionCollection.h"

#include "bxcommon.h"

namespace opt {

  static std::string bam; // the bam to analyze
  static bool verbose = false; 
  static std::string tag = "BX";
}

static const char* shortopts = "hvt:";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "verbose",                 no_argument, NULL, 'v' },
  { "tag",                     required_argument, NULL, 't' },
  { NULL, 0, NULL, 0 }
};

static const char *MOL_USAGE_MESSAGE =
"Usage: bxtools mol <BAM/SAM/CRAM> > mol.bed\n"
"Description: Return span of molecules from 10X data (using MI tag)\n"
"\n"
"  General options\n"
"  -v, --verbose         Set verbose output\n"
"  -t, --tag             Use a different tag other than MI\n"
"\n";

class BXMol {
  
public:
  
  BXMol() {}

  int min = INT_MAX;
  int max = -1;
  int chr = -1;

  std::unordered_set<std::string> bx; // BX tag
  //std::string bx; // BX tag
  std::string mi; // MI (or -t) tag
  std::string tmpbx; // tmp to be overwritted to hold new bx
  
  int nr = 0; // num reads
  std::string chr_string;

  bool add(const SeqLib::BamRecord& r, const SeqLib::BamHeader h) {
    
    if (chr > 0 && chr != r.ChrID()) {
      std::cerr << "Warning: " << opt::tag << " "  << mi << " spans multiple chromosomes" << std::endl;
      return false;
    }

    
    r.GetTag(opt::tag, mi);

    ++nr;
    
    // get the BX tag
    r.GetTag("BX", tmpbx);
    bx.insert(tmpbx);

    // set the position
    chr = r.ChrID();
    min = std::min(r.Position(), min);
    max = std::max(r.PositionEnd(), max);
    if (chr_string.empty())
      chr_string = h.IDtoName(chr);
    
    return true;

  }
  
  friend std::ostream& operator<<(std::ostream& out, const BXMol& b) {
    std::stringstream ss;
    for (auto& i : b.bx)
      ss << i << ",";
    std::string bxstring = ss.str();
    if (!bxstring.empty())
      bxstring.pop_back(); // remove last comma
    out << b.chr_string << "\t" << b.min << "\t" 
	<< b.max << "\t" << b.mi << "\t" << bxstring << "\t" 
	<< b.nr;
    return out;
  }

};

static void parseOptions(int argc, char** argv);

void runMol(int argc, char** argv) {
  
  parseOptions(argc, argv);

  SeqLib::BamReader reader;
  BXOPEN(reader, opt::bam);
  SeqLib::BamHeader hdr = reader.Header();

  std::unordered_map<std::string, BXMol> molmap;

  SeqLib::BamRecord r;
  size_t count = 0; 
  std::string mi;
  //int32_t mi;
  while (reader.GetNextRecord(r)) {
    BXLOOPCHECK(r, molmap.size(), opt::tag);
    if (r.MappedFlag() && r.GetTag(opt::tag, mi)) 
      molmap[mi].add(r, hdr);
  }  
  // print them out as a BED
  for (const auto& b : molmap)
    std::cout << b.second << std::endl;
}

static void parseOptions(int argc, char** argv) {

  bool die = false;
  bool help = false;

  if (argc < 2) 
    die = true;
  else
    opt::bam = std::string(argv[1]);

  std::stringstream ss;

  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v': opt::verbose = true; break;
    case 'h': help = true; break;
    case 't': arg >> opt::tag; break;
    }
  }

  if (die || help) {
    std::cerr << "\n" << MOL_USAGE_MESSAGE;
    die ? exit(EXIT_FAILURE) : exit(EXIT_SUCCESS);
  }
  
}

