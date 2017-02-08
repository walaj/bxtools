#ifndef BXTOOLS_STATS_H__
#define BXTOOLS_STATS_H__

#include <cstdint>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>

void runStat(int argc, char** argv);

struct BXStat {

  std::string bx; // label
  size_t count;   // number of reads
  std::vector<int> isize; // insert size
  std::vector<int> mapq;  // mapping quality
  
  friend std::ostream& operator<<(std::ostream& out, const BXStat& b);

};


#endif
