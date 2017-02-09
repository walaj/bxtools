#include "bxconvert.h"

#include <getopt.h>
#include <iostream>
#include <sstream>

#include "SeqLib/BamReader.h"
#include "SeqLib/BamWriter.h"
#include "SeqLib/GenomicRegionCollection.h"

#include "bxcommon.h"
#include <cassert>

namespace opt {

  static std::string bam; // the bam to analyze
  static std::string analysis_id = "foo"; // unique prefix for new BAM output
  static bool verbose = false;
  static int width = 1000;
  static std::string bed; // optional bed file                                  
}

static const char* shortopts = "hvw:a:b:";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "bed",                     required_argument, NULL, 'b' },
  { "pad",                     required_argument, NULL, 'p' },
  { "width",                   required_argument, NULL, 'w' },
  { "analysis-id",                 required_argument, NULL, 'a' },
  { NULL, 0, NULL, 0 }
};


  void runConvert(int argc, char** argv) {

    parseConvertOptions(argc, argv);

    SeqLib::BamReader reader;
    BXOPEN(reader, opt::bam);
    SeqLib::BamHeader hdr = reader.Header();
    
    SeqLib::BamRecord r;
    SeqLib::BamWriter w;
    size_t count = 0;
    std::unordered_map<std::string, size_t> bxtags;
    std::stringstream ss;
    const std::string empty_tag = "Empty";
    if (opt::bam.compare("-") == 0){
      std::cerr << "Cant accept standard input as file" << std::endl;
      exit(EXIT_FAILURE);
    }
    // Loop through file once to grab all BX tags and store in string to generate header
    ss << "@HD" << "\t" << "VN:1.4" << "  " << "GO:none SO:unsorted" << std::endl;  
    while (reader.GetNextRecord(r)){
      std::string bx = r.GetZTag("BX");
      if (bx.empty()){
	bx = empty_tag;
      }
      assert(!bx.empty());
      if (!bxtags.count(bx)) {
        bxtags.insert(std::pair<std::string, size_t>(bx, count));
	++count;      
	ss << "@SQ" << "\t" << "SN:" << bx << "\t" << "LN:500000000" << std::endl;
      }    
      
    }
    //write new header based on string generated from BX tags
    SeqLib::BamHeader bxbamheader (ss.str());
    //std::string newbname = opt::analysis_id + "_BXsorted" + ".bam";
    w.Open("/broad/hptmp/tkamath/foobar.bam");
    w.SetHeader(bxbamheader);
    w.WriteHeader();
    
    //Loop through the BAM file again
    reader.Close();
    SeqLib::BamReader reader2;
    BXOPEN(reader2, opt::bam);
    SeqLib::BamRecord rr;
    while (reader2.GetNextRecord(rr)){
      std::string bx = rr.GetZTag("BX");
      if (bx.empty()){
	bx = empty_tag;
      }

        std::string chr = hdr.IDtoName(rr.ChrID());
        rr.SetChrID(bxtags[bx]);
        rr.AddZTag("CR", chr); 
        w.WriteRecord(rr);
    }
    w.Close();
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
      case 'a': arg >> opt::analysis_id; break;
      case 'b': arg >> opt::bed; break;
      }
    }

}
