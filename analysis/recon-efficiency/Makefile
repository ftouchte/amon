#********************************************************
# MAKEFILE for BUILDING EXAMPLES FOR HIPO4 LIBRARY
# AUTHOR: GAVALIAN DATE: 10/24/2018
#
# Version modifiée par Felix Touchte Codjo pour son 
# usage personnel. Merci GAVALIAN ! 
# Date : 12 avril 2024
#********************************************************

PATH2HIPO := /home/touchte-codjo/framework/hipo

HIPOCFLAGS  := -I$(PATH2HIPO)/hipo4 -I$(PATH2HIPO)/hipo4/chart   
HIPOLIBS    := -L$(PATH2HIPO)/lib -lhipo4 
									
LZ4LIBS     := -L$(PATH2HIPO)/lz4/lib -llz4
LZ4INCLUDES := -I$(PATH2HIPO)/lz4/lib


# ROOT libraries 
ROOTLIBS = $(shell root-config --libs)
# ROOT include flags
ROOTCFLAGS = $(shell root-config --cflags)

GTKLIBS = $(shell pkg-config --libs gtkmm-4.0)
GTKFLAGS = $(shell pkg-config --cflags gtkmm-4.0)

CAIROLIBS = $(shell pkg-config --libs cairomm-1.16)
CAIROFLAGS = $(shell pkg-config --cflags cairomm-1.16)


CXX       := g++
CXXFLAGS  += -Wall -fPIC -std=c++17
LD        := g++
LDFLAGS   :=


#all:  showFile histo plot benchmark simu
all: deep_process occupancy elastics

occupancy: occupancy.o
	$(CXX) -o occupancy.exe $^ $(HIPOLIBS) $(LZ4LIBS) $(ROOTLIBS)

deep_process: deep_process.o futils.o
	$(CXX) -o deep_process.exe $^ $(HIPOLIBS) $(LZ4LIBS) $(ROOTLIBS)

elastics: elastics.o futils.o
	$(CXX) -o elastics.exe $^ $(HIPOLIBS) $(LZ4LIBS) $(ROOTLIBS)


# $< représente la première de la cible, i.e histo.o
# $^ représente la liste complète des dépendances

clean:
	@echo 'Removing all build files'
	@rm -rf *.o *~ *.exe example*hipo ./output/*.pdf

%.o: %.cpp
	$(CXX) -c $< -O2 $(CXXFLAGS) $(HIPOCFLAGS) $(LZ4INCLUDES) $(ROOTCFLAGS) $(GTKFLAGS) $(CAIROFLAGS) 
	

