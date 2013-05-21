TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += main.cpp

QMAKE_CXXFLAGS += -std=c++11
QMAKE_LIBDIR += /usr/local/lib

LIBS += -lboost_program_options-mt

HEADERS += \
    solutions.h \
    pre_integration.h \
    integration.h
