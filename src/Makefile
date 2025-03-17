#********************************************************
# MAKEFILE for BUILDING EXAMPLES FOR HIPO4 LIBRARY
# AUTHOR: GAVALIAN DATE: 10/24/2018
#
# Version modifiée par Felix Touchte Codjo pour son 
# usage personnel. Merci GAVALIAN ! 
# Date : 12 avril 2024
#********************************************************

# ATTENTION

# Les "paths" sont en chemin relatifs (rq: les ".."). Prière de localiser le dossier *hipo/hipo4...
# Quelques options de compilations
# 	-I<path> : Ajoute le répertoire <path> au chemin de recherche des fichiers d'en-tête.
# 	-L<path> : Ajoute le répertoire <path> au chemin de recherche des bibliothèques.

#PATH2HIPO := /homeijclab/touchte-codjo/hipo
#PATH2HIPO := /home/ftouchte/hipo
PATH2HIPO := /home/touchte-codjo/framework/hipo

HIPOCFLAGS  := -I$(PATH2HIPO)/hipo4 -I$(PATH2HIPO)/hipo4/chart   
HIPOLIBS    := -L$(PATH2HIPO)/lib -lhipo4 
									
LZ4LIBS     := -L$(PATH2HIPO)/lz4/lib -llz4
LZ4INCLUDES := -I$(PATH2HIPO)/lz4/lib

FELIXFLAGS := -I/homeijclab/touchte-codjo/Bureau/alert/cpp/ahdc

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
all: gui  


histAhdcAdc: histAhdcAdc.o AhdcExtractor.o 
	$(CXX) -o histAhdcAdc.exe $^ $(HIPOLIBS) $(LZ4LIBS) $(ROOTLIBS) 

test_fH1D: test_fH1D.o fAxis.o fH1D.o fCanvas.o fH2D.o
	$(CXX) -o test_fH1D.exe $^ $(CAIROLIBS)

test_fAxis: test_fAxis.o fAxis.o 
	$(CXX) -o test_fAxis.exe $^

gui: gui.o AhdcExtractor.o AhdcDetector.o Point3D.o fAxis.o fCanvas.o fH1D.o fColorPalette.o
	$(CXX) -o gui.exe $^ $(HIPOLIBS) $(LZ4LIBS) $(ROOTLIBS) $(GTKLIBS)



# $< représente la première de la cible, i.e histo.o
# $^ représente la liste complète des dépendances

clean:
	@echo 'Removing all build files'
	@rm -rf *.o *~ *.exe example*hipo *.pdf

%.o: %.cpp
	$(CXX) -c $< -O2 $(CXXFLAGS) $(HIPOCFLAGS) $(LZ4INCLUDES) $(ROOTCFLAGS) $(GTKFLAGS) 
	
#	$(CXX) -c $< -O2 $(CXXFLAGS) $(HIPOCFLAGS) $(LZ4INCLUDES) $(ROOTCFLAGS) $(RAPID_CHECK_FLAGS) $(FELIXFLAGS)

