#include "bxtile.h"

#include <getopt.h>
#include <iostream>
#include <sstream>

#include "SeqLib/BamReader.h"
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

static const char *TILE_USAGE_MESSAGE =
"Usage: bxtools tile <BAM/SAM/CRAM> > tiles.bed\n"
"Description: Gather BX counts on tiled ranges\n"
"\n"
"  General options\n"
"  -v, --verbose         Set verbose output\n"
"  -w, --width           Width of the tile [1000]\n"
"  -O, --overlap         Overlap of the tiles [0]\n"
"  -b, --bed             Rather than tile genome, input BED with regions\n"
"\n";

class BXRegion : public SeqLib::GenomicRegion {
  
public:
  
  BXRegion() : GenomicRegion() {}

  BXRegion(const std::string c, const std::string p1, const std::string p2, 
	   const SeqLib::BamHeader& h) : GenomicRegion(c, p1, p2, h) {}

  std::unordered_map<std::string, size_t> counts;

  std::string ToBEDString(const SeqLib::BamHeader& h) const {
    std::string out = h.IDtoName(chr) + "\t" + std::to_string(pos1) + 
      "\t" + std::to_string(pos2);
    if (counts.size())
      out += "\t";
    for (const auto& b : counts)
      out +=  b.first + "_" + std::to_string(b.second) + ",";
    if (counts.size())
      out.pop_back(); // erase last comma
    return out;
  }
};

void runTile(int argc, char** argv) {
  
  parseTileOptions(argc, argv);

  SeqLib::BamReader reader;
  BXOPEN(reader, opt::bam);
  SeqLib::BamHeader hdr = reader.Header();

  SeqLib::GenomicRegionCollection<BXRegion> * tiles = nullptr;
  if (!opt::bed.empty()) {
    tiles = new SeqLib::GenomicRegionCollection<BXRegion>();
    tiles->ReadBED(opt::bed, hdr);
    tiles->CreateTreeMap();
  } else {
    // tile it
    std::cerr << "...creating tiles with width " << 
      SeqLib::AddCommas(opt::width) << " and overlap " << SeqLib::AddCommas(opt::overlap) << std::endl;
    tiles = new SeqLib::GenomicRegionCollection<BXRegion>(opt::width, opt::overlap, hdr.GetHeaderSequenceVector());
    std::cerr << "...created " << SeqLib::AddCommas(tiles->size()) << " tiles" << std::endl;
    std::cerr << "...sorting and creating interval tree" << std::endl;
    tiles->CreateTreeMap();
  }

  std::cerr << "...reading input" << std::endl;
  SeqLib::BamRecord r;
  size_t count = 0; 
  size_t bxcount = 0;
  while (reader.GetNextRecord(r)) {
    BXLOOPCHECK(r, bxcount);

    if (r.MappedFlag()) {
      std::vector<int> bins = tiles->FindOverlappedIntervals(r.AsGenomicRegion(), true);
      for (const auto& b : bins) 
	++(*tiles)[b].counts[bx];
      ++bxcount;
    }
      
  }

  for (const auto& b : *tiles)
    std::cout << b.ToBEDString(hdr) << std::endl;

  if (tiles)
    delete tiles;
  
}

void parseTileOptions(int argc, char** argv) {

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
    case 'w': arg >> opt::width; break;
    case 'O': arg >> opt::overlap; break;
    case 'b': arg >> opt::bed; break;
    }
  }

  if (die || help) {
    std::cerr << "\n" << TILE_USAGE_MESSAGE;
    die ? exit(EXIT_FAILURE) : exit(EXIT_SUCCESS);
  }
  
}

