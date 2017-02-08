#ifndef BXTOOLS_BXCOMMON_H__
#define BXTOOLS_BXCOMMON_H__

#define BXOPEN(reader, bam)			\
  if (!reader.Open(bam)) {				     \
    std::cerr << "Failed to open bam: " << bam << std::endl; \
    exit(EXIT_FAILURE); \
  }			\

#define BXLOOPCHECK(r, found, tag)						\
    ++count;								\
    if (count == 100000 && !(found))					\
      std::cerr << "****1e5 reads in and haven't hit " << tag << " tag yet****" << std::endl; \
    if (count == 1000000 && !(found))					\
      std::cerr << "****1e6 reads in and haven't hit " << tag << " tag yet****" << std::endl; \
    if (count % 1000000 == 0 && opt::verbose)					\
      std::cerr << "...at read " << SeqLib::AddCommas(count) << " at pos " << r.Brief() << std::endl;
    
#endif
