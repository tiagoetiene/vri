TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += main.cpp \
    GageAdaptor.cpp

QMAKE_CXXFLAGS += -std=c++11
QMAKE_LIBDIR += /usr/local/lib

LIBS += -lboost_program_options-mt -lteem

HEADERS += \
    solutions.h \
    pre_integration.h \
    integration.h \
    io.h \
    GageAdaptor.h
