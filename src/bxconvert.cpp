#include "bxconvert.h"

#include <getopt.h>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <cassert>

#include "SeqLib/BamReader.h"
#include "SeqLib/BamWriter.h"
#include "SeqLib/GenomicRegionCollection.h"

#include "bxcommon.h"

static const char *CONVERT_USAGE_MESSAGE =
"Usage: bxtools convert <BAM/SAM/CRAM> > converted.bam\n"
"Description: Convert a BAM to a BX sorted BAM by switching BX and chromosome\n"
"\n"
"  General options\n"
"  -v, --verbose         Set verbose output\n"
"  -k, --keep-tags       Add chromosome tag (CR) and position (PS) tag, and keep other tags. Default: delete all tags\n"
"  -t, --tag             Tag to flip for chromosome. Default: BX\n"
"\n";

namespace opt {
  static bool verbose = false;
  static std::string bam;
  static bool keeptags = false;
  static std::string tag = "BX";
}

static const char* shortopts = "hvkt:";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "verbose",                 no_argument, NULL, 'v' },
  { "keep-tags",               no_argument, NULL, 'k' },
  { "tag",                     required_argument, NULL, 't' },
  { NULL, 0, NULL, 0 }
};

static void read_bx(std::string& bx, const SeqLib::BamRecord& r);
static const std::string empty_tag = "Empty";

void runConvert(int argc, char** argv) {

    parseOptions(argc, argv);

    SeqLib::BamReader reader;
    BXOPEN(reader, opt::bam);
    SeqLib::BamHeader hdr = reader.Header();
    
    SeqLib::BamRecord r;
    SeqLib::BamWriter w;
    size_t count = 0, unique_bx = 0;
    std::string bx;
    std::unordered_map<std::string, size_t> bxtags;
    std::stringstream ss;

    if (opt::bam.compare("-") == 0){
      std::cerr << "Cant accept standard input as file" << std::endl;
      exit(EXIT_FAILURE);
    }

    if (opt::verbose)
      std::cerr << "...starting first pass to tally unique " << opt::tag << " tags" << std::endl;

    // Loop through file once to grab all BX tags and store in string to generate header
    ss << "@HD" << "\t" << "VN:1.4" << "  " << "GO:none\tSO:unsorted" << std::endl;  
    while (reader.GetNextRecord(r)){
      read_bx(bx, r);

      BXLOOPCHECK(r, unique_bx > 1, opt::tag)
      if (!bxtags.count(bx)) {
        bxtags.insert(std::pair<std::string, size_t>(bx, unique_bx));
	++unique_bx;      
	ss << "@SQ" << "\t" << "SN:" << bx << "\t" << "LN:1" << std::endl;
      }    
      
    }
    //write new header based on string generated from BX tags
    SeqLib::BamHeader bxbamheader (ss.str());

    if (opt::verbose) 
      std::cerr << "Found " << unique_bx << " unique barcodes" << std::endl;

    w.Open("-");
    w.SetHeader(bxbamheader);
    w.WriteHeader();
    
    //Loop through the BAM file again
    reader.Close();
    SeqLib::BamReader reader2;
    BXOPEN(reader2, opt::bam);
    
    if (opt::verbose)
      std::cerr << "...starting second pass to flip chr and BX" << std::endl;

    while (reader2.GetNextRecord(r)) {

      if (opt::keeptags) {
	r.AddZTag("CR", r.ChrID() >= 0 ? hdr.IDtoName(r.ChrID()) : "*");
	r.AddIntTag("POS", r.Position());
      }

      read_bx(bx, r); // read the BX tag. Set default if not present
      BXLOOPCHECK(r, true, opt::tag) // read and check we have a BX
      r.SetChrID(bxtags[bx]);
      r.SetChrIDMate(-1);
      r.SetPosition(0);
      if (!opt::keeptags)
	r.RemoveAllTags();

      w.WriteRecord(r);
    }
    w.Close();
  }


void parseOptions(int argc, char** argv) {

  bool die = false;
  bool help = false;

  if (argc < 2)
    die = true;
  else
    opt::bam = std::string(argv[1]);

  std::stringstream ss;

  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)	\
    {
      std::istringstream arg(optarg != NULL ? optarg : "");
      switch (c) {
      case 'v': opt::verbose = true; break;
      case 'h': help = true; break;
      case 'k': opt::keeptags = true; break;
      case 't': arg >> opt::tag;  break;
      }
    }

  if (die || help) {
    std::cerr << "\n" << CONVERT_USAGE_MESSAGE;
    die ? exit(EXIT_FAILURE) : exit(EXIT_SUCCESS);
  }


}

static void read_bx(std::string& bx, const SeqLib::BamRecord& r) {
  if (!r.GetZTag(opt::tag, bx))
    bx = empty_tag;
  std::replace(bx.begin(), bx.end(), '-', '_');
  assert(!bx.empty());
}
