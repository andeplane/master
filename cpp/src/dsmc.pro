TEMPLATE = app
CONFIG += console
CONFIG -= qt
LIBS   += -larmadillo -openmp
INCLUDEPATH += -openmp
QMAKE_CXXFLAGS += -openmp
QMAKE_LFLAGS += -openmp

SOURCES += main.cpp \
    system.cpp \
    cell.cpp \
    molecule.cpp \
    sorter.cpp \
    statisticssampler.cpp \
    wall.cpp \
    Random.cpp \
    grid.cpp

HEADERS += \
    system.h \
    cell.h \
    molecule.h \
    sorter.h \
    statisticssampler.h \
    wall.h \
    Random.h \
    grid.h

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../CiRIOBase/release/ -lCiRIOBase
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../CiRIOBase/debug/ -lCiRIOBase
else:symbian: LIBS += -lCiRIOBase
else:unix: LIBS += -L$$PWD/../CiRIOBase/ -lCiRIOBase

INCLUDEPATH += $$PWD/../CiRIOBase/source
DEPENDPATH += $$PWD/../CiRIOBase/source

win32:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../CiRIOBase/release/CiRIOBase.lib
else:win32:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../CiRIOBase/debug/CiRIOBase.lib
else:unix:!symbian: PRE_TARGETDEPS += $$PWD/../CiRIOBase/libCiRIOBase.a

OTHER_FILES += \
    ../dsmc.ini
