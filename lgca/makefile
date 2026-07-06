LPATH=
IPATH= -I$(HOME)/Documents/lgca/rng/
INCL= $(HOME)/Documents/lgca/rng/WELL1024a.c
LINK= -lm 
CXXFLAGS=-O2

CIMGLINK= -lpthread -lX11
CIMGPATH= -L/usr/X11R6/lib

all:
	c++ -std=c++11 lgca_noib.cpp $(INCL) $(LPATH) $(LINK) $(CXXFLAGS) -o nem2b_print
# 	c++ -std=c++11 base_rocket.cpp $(INCL) $(LPATH) $(LINK) $(CXXFLAGS) -o base_woundout
# 	c++ -std=c++11 lgca_pdediscr.cpp $(INCL) $(LPATH) $(LINK) $(CXXFLAGS) -o phasesp_sweep
#	c++ -std=c++11 base_gsqhex.cpp $(INCL) $(LPATH) $(LINK) $(CXXFLAGS) -o persnem_sweep
all_time: 
	c++ base-time.cpp $(INCL) $(LPATH) $(LINK) $(CXXFLAGS) -o base-time
all_cimg:
	c++ -std=c++11 lgca_clean-1.cpp -DCIMG=1 $(CIMGLINK) $(CIMGPATH) $(INCL) $(LPATH) $(LINK) $(CXXFLAGS) -g -o lgca_test
#	c++ -std=c++11 lgca_noib.cpp -DCIMG=1 $(CIMGLINK) $(CIMGPATH) $(INCL) $(LPATH) $(LINK) $(CXXFLAGS) -g -o lgca_poln
#	c++ -std=c++11 lgca_nenmcorr.cpp -DCIMG=1 $(CIMGLINK) $(CIMGPATH) $(INCL) $(LPATH) $(LINK) $(CXXFLAGS) -g -o lgca_ncorrimgs
# 	c++ -std=c++11 base_rocket.cpp -DCIMG=1 $(CIMGLINK) $(CIMGPATH) $(INCL) $(LPATH) $(LINK) $(CXXFLAGS) -g -o base_wound
#	c++ -std=c++11 base_gsqhex.cpp -DCIMG=1 $(CIMGLINK) $(CIMGPATH) $(INCL) $(LPATH) $(LINK) $(CXXFLAGS) -g -o base_imgs
#	c++ -std=c++11 base_boson.cpp -DCIMG=1 $(CIMGLINK) $(CIMGPATH) $(INCL) $(LPATH) $(LINK) $(CXXFLAGS) -g -o base_test
# 	c++ -std=c++11 base_imgcollect.cpp -DCIMG=1 $(CIMGLINK) $(CIMGPATH) $(INCL) $(LPATH) $(LINK) $(CXXFLAGS) -g -o base_gnemat
all_half:
	c++ base-half.cpp $(INCL) $(LPATH) $(LINK) $(CXXFLAGS) -o base-half

all_para:
	mpic++ base.cpp -DPARA=1 $(INCL) $(LPATH) $(LINK) $(CXXFLAGS) -o base
debug:
	c++ -std=c++11 lgca_nempers.cpp $(INCL) $(LPATH) $(LINK) -g -o base_testnoib
#	c++ base_dbug.cpp $(INCL) $(LPATH) $(LINK) -g -o base_bug
# 	c++ base_dbug.cpp $(INCL) $(LPATH) $(LINK) $(CXXFLAGS) -ggdb -o base_bug
prof:
	c++ base.cpp $(INCL) $(LPATH) $(LINK) -pg $(CXXFLAGS) -o base
all_david:
	c++ -std=c++11 base_imgdavid.cpp -DCIMG=1 $(CIMGLINK) $(CIMGPATH) $(INCL) $(LPATH) $(LINK) $(CXXFLAGS) -o base_david
debug_david:
	c++ -std=c++11 base_imgdavid.cpp $(INCL) $(LPATH) $(LINK) $(CXXFLAGS) -ggdb -o base_david
