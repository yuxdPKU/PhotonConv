##############################################
# please add new classes in alphabetical order

AUTOMAKE_OPTIONS = foreign

# List of shared libraries to produce
lib_LTLIBRARIES = \
  libtrack_to_calo.la

#
libtrack_to_calo_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib64 \
  -L$(OFFLINE_MAIN)/lib

libtrack_to_calo_la_LIBADD = \
  -lActsCore \
  -lActsPluginTGeo \
  -lActsExamplesCommon \
  -lActsExamplesDetectorTGeo \
  -lActsExamplesFramework \
  -lSubsysReco \
  -ltrackbase_historic_io \
  -ltpc \
  -ltrack_io \
  -ltrack_reco \
  -lcalo_io \
  -lglobalvertex_io \
  -lg4detectors_io \
  -lphparameter \
  -ldecayfinder_io \
  -lphhepmc \
  -lkfparticle_sphenix_io \
  -lkfparticle_sphenix

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include \
  -DHomogeneousField

pkginclude_HEADERS = \
  SeedContainerMaker.h \
  TrackCaloMatch.h \
  TrackToCalo.h \
  CaloOnly.h \
  TrackOnly.h

libtrack_to_calo_la_SOURCES = \
  SeedContainerMaker.cc \
  TrackCaloMatch.cc \
  TrackToCalo.cc \
  CaloOnly.cc \
  TrackOnly.cc

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals_track_to_calo

BUILT_SOURCES  = testexternals.cc

testexternals_track_to_calo_SOURCES = testexternals.cc
testexternals_track_to_calo_LDADD = libtrack_to_calo.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

##############################################
# please add new classes in alphabetical order

clean-local:
	rm -f $(BUILT_SOURCES)
