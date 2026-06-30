# Higgstools-in-2HDMC
edit to 2HDMC code to include updated limits from HiggsTools

***This might not be optimised, and i cannot promise zero bugs.***

2HDMC code can be found at [2HDMC](https://2hdmc.hepforge.org/)

HiggsBounds, HiggsTools, and the datasets for HiggsBounds and HiggsSignals can be found at [HiggsBounds](https://gitlab.com/higgsbounds)


## step 1 --- Update Makefile

edit the Makefile to include the lines instead of the import of HiggsBounds and HiggsSignals
```
# To use HiggsBounds/HiggsTools for Higgs constraints, set both of the
# following path variables to the corresponding build directories.
# Requires HiggsBounds>=5.7.0 and HiggsSignals>=2.4.0
 HiggsTools_DIR="path-to"/higgstools/build
 HiggsBounds_DIR="path-to"/higgsbounds/build

ifdef HiggsBounds_DIR
ifdef HiggsTools_DIR
CFLAGS+=-DHiggsBounds
LDFLAGS+=-L$(HiggsBounds_DIR)/lib -lHB -lgfortran
LDFLAGS+=-L$(HiggsTools_DIR)/src -lHiggsTools
INCLUDE+=-I$(HiggsTools_DIR)/../include
INCLUDE+=-I$(HiggsBounds_DIR)/../include
SOURCES+=HBHS.cpp
endif
endif
```
Note HiggsBounds is still used to find cross sections for ${H^\pm}$ production at the LHC.

## step 2 --- Include new HBHS file
use the new `HBHS.h` and `HBHS.cpp` files in the `./src` directory of 2HDMC.

## step 3 --- Update paths to HiggsBounds and HiggsSignals datasets in `HBHS.h`
edit lines `13-14` of `HBHS.h` to include the correct path to the HiggsBounds / HiggsSignals datasets.
