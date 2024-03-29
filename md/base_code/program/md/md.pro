TEMPLATE = app
CONFIG += console
CONFIG -= qt
TARGET = main

release {

}

SOURCES += \
    system.cpp \
    statisticssampler.cpp \
    md.cpp \
    random.cpp \
    thermostat.cpp \
    cutil.cpp \
    unitconverter.cpp \
    settings.cpp \
    mdio.cpp \
    mdtimer.cpp

HEADERS += \
    system.h \
    statisticssampler.h \
    random.h \
    thermostat.h \
    cutil.h \
    cinifile.h \
    unitconverter.h \
    settings.h \
    mdio.h \
    mdtimer.h \
    atom_types.h \
    potential_lennard_jones.h

mac {
    QMAKE_CXX = icpc
    CONFIG -= app_bundle
    LIBS   +=
    INCLUDEPATH +=
    QMAKE_CXXFLAGS +=
    QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS
    QMAKE_CXXFLAGS_DEBUG = $$QMAKE_CXXFLAGS
}

unix:!mac {

    QMAKE_CXX = icpc
    LIBS   +=
    INCLUDEPATH +=
    QMAKE_CXXFLAGS +=
    QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS
    QMAKE_CXXFLAGS_DEBUG = $$QMAKE_CXXFLAGS
}

# MPI Settings
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX

QMAKE_CFLAGS = $$system(mpicc --showme:compile)
QMAKE_LFLAGS = $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS
QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE += -O3

QMAKE_LFLAGS += -ipo -no-prec-div -falign-functions=16 -xCORE-AVX-I
QMAKE_CXXFLAGS_RELEASE += -ipo -no-prec-div -falign-functions=16 -xCORE-AVX-I

#QMAKE_CXXFLAGS_RELEASE += -ipo -prof-use -no-prec-div -falign-functions=16 -xCORE-AVX-I
#QMAKE_CXXFLAGS_RELEASE += -ipo -prof-gen
#QMAKE_LFLAGS += -falign-functions=16 -xCORE-AVX-I -ipo -prof-gen

QMAKE_LFLAGS -= -lm
QMAKE_LFLAGS -= -O2
QMAKE_LFLAGS += -O3
QMAKE_CFLAGS_RELEASE -= -fPIE
QMAKE_CXXFLAGS_RELEASE -= -fPIE
