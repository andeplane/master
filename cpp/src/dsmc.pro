TEMPLATE = app
CONFIG += console
CONFIG -= qt
LIBS   += -larmadillo

SOURCES += main.cpp \
    system.cpp \
    cell.cpp \
    collisionobject.cpp \
    molecule.cpp \
    sorter.cpp \
    statisticssampler.cpp \
    wall.cpp \
    Random.cpp

HEADERS += \
    system.h \
    cell.h \
    collisionobject.h \
    molecule.h \
    sorter.h \
    statisticssampler.h \
    wall.h \
    Random.h

