TEMPLATE = app
DESTDIR = bin/

CONFIG = debug_and_release release
QMAKE_CXXFLAGS += -fopenmp
QMAKE_LFLAGS += -fopenmp
QT =

FREEDLIBS = $(HOME)/src/cpp/freedlibs/

INCLUDEPATH = \
  $$FREEDLIBS/libexception/src \
  $$FREEDLIBS/libaline/src \
  $$FREEDLIBS/libgeneral/src \
  $$FREEDLIBS/libshape/src \
  $$FREEDLIBS/libimage/src \
  $$FREEDLIBS/libfreed/src

LIBS = \
  -L$$FREEDLIBS/libexception/src \
  -L$$FREEDLIBS/libaline/src  \
  -L$$FREEDLIBS/libgeneral/src \
  -L$$FREEDLIBS/libshape/src  \
  -L$$FREEDLIBS/libimage/src  \
  -L$$FREEDLIBS/libfreed/src \
#  -lfreed -limage -lshape -lgeneral -laline -lexception -ltiff
  -limage -lshape -lgeneral -laline -lexception -ltiff

SOURCES += \
    src/main.cpp \
    src/findmicrotubules.cpp \
    src/findparticles.cpp \
    src/findparticlesmanually.cpp \
    src/analysingparticles.cpp \
    src/getobjectinterdistances.cpp
