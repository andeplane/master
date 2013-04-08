TEMPLATE = app
CONFIG += console
CONFIG -= qt

release {
    DEFINES +=
}

DEFINES +=

SOURCES += main.cpp \
    system.cpp \
    cell.cpp \
    molecule.cpp \
    statisticssampler.cpp \
    random.cpp \
    grid.cpp \
    cutil.cpp \
    system.inc.cpp \
    unitconverter.cpp \
    image.cpp \
    settings.cpp \
    dsmc_io.cpp \
    threadcontrol.cpp \
    dsmctimer.cpp

HEADERS += \
    system.h \
    cell.h \
    molecule.h \
    statisticssampler.h \
    random.h \
    grid.h \
    cutil.h \
    cinifile.h \
    unitconverter.h \
    image.h \
    defines.h \
    settings.h \
    dsmc_io.h \
    threadcontrol.h \
    dsmctimer.h

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
    QMAKE_CXX = mpic++
}

QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = mpicc
QMAKE_CFLAGS = $$system(mpicc --showme:compile)
QMAKE_LFLAGS = $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS
