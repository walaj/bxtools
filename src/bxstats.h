#ifndef BXTOOLS_STATS_H__
#define BXTOOLS_STATS_H__

#include <cstdint>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>

void parseStatOptions(int argc, char** argv);
void runStat(int argc, char** argv);

//http://stackoverflow.com/questions/2114797/compute-median-of-values-stored-in-vector-c
static double CalcMHWScore(std::vector<int> scores);

struct BXStat {

  std::string bx; // label
  size_t count;   // number of reads
  std::vector<int> isize; // insert size
  std::vector<int> mapq;  // mapping quality
  
  friend std::ostream& operator<<(std::ostream& out, const BXStat& b) {
    double isize_med = -1;
    double mapq_med = -1;
    if (b.isize.size())
      isize_med = CalcMHWScore(b.isize);
    if (b.mapq.size())
      mapq_med = CalcMHWScore(b.mapq);
    out << b.bx << "\t" << b.count << "\t" << isize_med << "\t" << mapq_med;
    return out;
  }

};


#endif
