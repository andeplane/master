TEMPLATE = app
CONFIG += console
CONFIG -= qt

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
    Image.cpp

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
    Image.h

OTHER_FILES += \
    ../dsmc.ini

mac {
    LIBS   += -larmadillo -openmp
    INCLUDEPATH +=
    QMAKE_CXXFLAGS += -openmp
    QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS
    QMAKE_CXXFLAGS_DEBUG = $$QMAKE_CXXFLAGS
}

unix:!mac {
    LIBS   += -larmadillo -fopenmp
    INCLUDEPATH +=
    QMAKE_CXXFLAGS += -fopenmp
    QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS
    QMAKE_CXXFLAGS_DEBUG = $$QMAKE_CXXFLAGS
}