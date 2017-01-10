bin_PROGRAMS = bxtools

bxtools_CPPFLAGS = \
     -I$(top_srcdir)/SeqLib \
     -I$(top_srcdir)/SeqLib/htslib -Wno-sign-compare

bxtools_LDADD = \
	$(top_builddir)/SeqLib/src/libseqlib.a \
	$(top_builddir)/SeqLib/htslib/libhts.a 

bxtools_SOURCES = bxtools.cpp bxsplit.cpp bxstats.cpp bxtile.cpp bxrelabel.cpp

