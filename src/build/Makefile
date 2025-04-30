CC = g++
CFLAGS = -Wall -Wextra -Wshadow -fno-stack-protector -static-libgcc -fopenmp -static-libstdc++ -std=c++20 

TARGETDIR = ../../bin/
TARGETNAME ?= SALI3D
TARGET = $(TARGETDIR)$(TARGETNAME)
SRCDIR = ../
LIBDIR = ../../lib
OBJDIR = ./

BOOSTDIR = C:/BOOST/boost_1_87_0
INCLUDES = -I$(LIBDIR)
INCLUDES += -I$(BOOSTDIR)

ifeq ($(TARGETNAME),SALI3D)
    # Source and object files for sali
    SRCS = $(SRCDIR)main_sali_3d.cpp $(LIBDIR)/crtbp.cpp 
    OBJS = $(OBJDIR)main_sali_3d.o $(OBJDIR)crtbp.o  
else ifeq ($(TARGETNAME),TRAJECTORY)
    # Source and object files for trajectory
    SRCS = $(SRCDIR)main_trajectory.cpp $(LIBDIR)/crtbp.cpp 
    OBJS = $(OBJDIR)main_trajectory.o $(OBJDIR)crtbp.o 
else ifeq ($(TARGETNAME),TDBG)
    # Source and object files for trajectory
    SRCS = $(SRCDIR)main_trajectory_debug.cpp $(LIBDIR)/crtbp.cpp 
    OBJS = $(OBJDIR)main_trajectory.o $(OBJDIR)crtbp.o  
else ifeq ($(TARGETNAME),SALI2D)
    # Source and object files for trajectory
    SRCS = $(SRCDIR)main_sali_2d.cpp $(LIBDIR)/crtbp.cpp 
    OBJS = $(OBJDIR)main_sali_2d.o $(OBJDIR)crtbp.o 
else ifeq ($(TARGETNAME), DE441)
    SRCS = $(SRCDIR)main_DE441_sali.cpp $(LIBDIR)/crtbp.cpp  $(LIBDIR)/chaosmap_generator.cpp 
    OBJS = $(OBJDIR)main_DE441_sali.o $(OBJDIR)crtbp.o $(OBJDIR)chaosmap_generator.o 
else ifeq ($(TARGETNAME), ZVC)
    # Source and object files for debug
    SRCS = $(SRCDIR)main_zvc.cpp
    OBJS = $(OBJDIR)main_zvc.o 
else ifeq ($(TARGETNAME), DEBUG1)
    # Source and object files for debug
    SRCS = $(SRCDIR)main_debug.cpp $(LIBDIR)/crtbp.cpp 
    OBJS = $(OBJDIR)main_debug.o $(OBJDIR)crtbp.o  
else ifeq ($(TARGETNAME), DEBUG2)
    # Source and object files for debug
    SRCS = $(SRCDIR)main_trajectorysearch.cpp $(LIBDIR)/crtbp.cpp
    OBJS = $(OBJDIR)main_trajectorysearch.o $(OBJDIR)crtbp.o  
endif

HEADERS = $(LIBDIR)/crtbp.hpp $(LIBDIR)/vector3d.hpp $(LIBDIR)/utils.hpp  $(LIBDIR)/chaosmap_generator.hpp

# Build rules
all: $(TARGET)

$(TARGET): $(OBJS) | $(TARGETDIR)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $^

$(OBJDIR)%.o: $(SRCDIR)%.cpp $(HEADERS) | $(OBJDIR)
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<

$(OBJDIR)%.o: $(LIBDIR)/%.cpp $(HEADERS) | $(OBJDIR)
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<

# Create directories
$(TARGETDIR) $(OBJDIR):
	mkdir -p $@

# Cleanup
clean:
	rm -f $(OBJDIR)*.o $(TARGET)

.PHONY: all clean