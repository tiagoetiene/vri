TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += main.cpp
QMAKE_CXXFLAGS += -std=c++11
QMAKE_LIBDIR += /usr/local/lib

HEADERS += \
    solutions.h
