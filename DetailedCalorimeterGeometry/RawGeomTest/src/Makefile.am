AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
    libcalogeomtest.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

pkginclude_HEADERS = \
  CaloGeomTest.h

libcalogeomtest_la_SOURCES = \
  CaloGeomTest.cc

# Perhaps a bit more than requested
libcalogeomtest_la_LIBADD = \
  -lcalo_io \
  -lphool \
  -lmbd_io \
  -lSubsysReco \
  -lcentrality_io \
  -lfun4all \
  -lg4dst \
  -lg4detectors_io \
  -lphg4hit


################################################
# linking tests

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.C
testexternals_LDADD = libcalogeomtest.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

clean-local:
	rm -f $(BUILT_SOURCES)
