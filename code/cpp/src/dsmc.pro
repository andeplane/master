TEMPLATE = app
CONFIG += console
CONFIG -= qt

release {
    DEFINES += ARMA_NO_DEBUG
}

SOURCES += main.cpp \
    system.cpp \
    cell.cpp \
    molecule.cpp \
    sorter.cpp \
    statisticssampler.cpp \
    Random.cpp \
    grid.cpp \
    Stdafx.cpp \
    CVector.cpp \
    CUtil.cpp \
    CMath.cpp \
    system.inc.cpp \
    CBitMap.cpp \
    unitconverter.cpp \
    Image.cpp \
    settings.cpp \
    dsmc_io.cpp \
    threadcontrol.cpp

HEADERS += \
    system.h \
    cell.h \
    molecule.h \
    sorter.h \
    statisticssampler.h \
    Random.h \
    grid.h \
    Stdafx.h \
    CVector.h \
    CUtil.h \
    CMatrix.h \
    CMath.h \
    CIniFile.h \
    CBitMap.h \
    unitconverter.h \
    Image.h \
    defines.h \
    settings.h \
    dsmc_io.h \
    threadcontrol.h

OTHER_FILES += \
    ../dsmc.ini

mac {
    CONFIG -= app_bundle
    LIBS   += -larmadillo -llapack -lblas
    INCLUDEPATH +=
    QMAKE_CXXFLAGS +=
    QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS
    QMAKE_CXXFLAGS_DEBUG = $$QMAKE_CXXFLAGS
    QMAKE_CXX = icc
}

unix:!mac {
    LIBS   += -larmadillo -llapack -lblas
    INCLUDEPATH +=
    QMAKE_CXXFLAGS +=
    QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS
    QMAKE_CXXFLAGS_DEBUG = $$QMAKE_CXXFLAGS
    QMAKE_CXX = g++-4.7
}


QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = mpicc
QMAKE_CFLAGS = $$system(mpicc --showme:compile)
QMAKE_LFLAGS = $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS
