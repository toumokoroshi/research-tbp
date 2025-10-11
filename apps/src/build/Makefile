# Compiler and flags
CC = g++
# Using -g for debugging, added -std=c++20
CFLAGS = -g -Wall -Wextra -Wshadow -fno-stack-protector -static-libgcc -fopenmp -static-libstdc++ -std=c++20 

# --- DIRECTORIES ---
# Use forward slashes for better compatibility
TARGETDIR = ../../bin/
SRCDIR = ../
LIBDIR = ../../lib
# OBJDIR is currently unused, but can be set to e.g. "obj/"
OBJDIR = 
# Use forward slashes for Boost path as well
BOOSTDIR = C:/Program_Languages/cpp_lib/boost-1_88_0

# --- TARGETS ---
TARGETNAME ?= SALI3D
# Correctly join paths
TARGET = $(TARGETDIR)$(TARGETNAME).exe

# --- INCLUDES ---
INCLUDES = -I$(LIBDIR)
INCLUDES += -I$(BOOSTDIR)

# --- SOURCE FILE SELECTION ---
# This section determines which files to compile based on the TARGETNAME variable
ifeq ($(TARGETNAME),SALI3D)
    # Source and object files for sali
    SRCS = $(SRCDIR)main_sali_3d.cpp $(LIBDIR)/crtbp.cpp 
    OBJS = $(OBJDIR)main_sali_3d.o $(OBJDIR)crtbp.o  
else ifeq ($(TARGETNAME),TRAJECTORY)
    # Source and object files for trajectory
    SRCS = $(SRCDIR)main_trajectory.cpp $(LIBDIR)/crtbp.cpp 
    OBJS = $(OBJDIR)main_trajectory.o $(OBJDIR)crtbp.o 
else ifeq ($(TARGETNAME),TDBG)
    # Source and object files for trajectory debug
    SRCS = $(SRCDIR)main_trajectory_debug.cpp $(LIBDIR)/crtbp.cpp 
    OBJS = $(OBJDIR)main_trajectory_debug.o $(OBJDIR)crtbp.o  
else ifeq ($(TARGETNAME),SALI2D)
    # Source and object files for sali 2d
    SRCS = $(SRCDIR)main_sali_2d.cpp $(LIBDIR)/crtbp.cpp 
    OBJS = $(OBJDIR)main_sali_2d.o $(OBJDIR)crtbp.o 
else ifeq ($(TARGETNAME), DE441)
    SRCS = $(SRCDIR)main_DE441_sali.cpp $(LIBDIR)/crtbp.cpp  $(LIBDIR)/chaosmap_generator.cpp 
    OBJS = $(OBJDIR)main_DE441_sali.o $(OBJDIR)crtbp.o $(OBJDIR)chaosmap_generator.o 
else ifeq ($(TARGETNAME), ZVC)
    # Source and object files for zvc
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

# --- HEADERS ---
# All header files your source files depend on
HEADERS = $(LIBDIR)/crtbp.hpp $(LIBDIR)/vector3d.hpp $(LIBDIR)/utils.hpp  $(LIBDIR)/chaosmap_generator.hpp

# ==============================================================================
# --- BUILD RULES ---
# ==============================================================================

# Default goal
all: $(TARGET)

# Rule to link the final executable
# This depends on all object files (.o) and requires the target directory to exist
$(TARGET): $(OBJS) | $(TARGETDIR)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $(OBJS)

# --- PATTERN RULES FOR COMPILING ---
# These rules tell make how to create .o files from .cpp files

# Rule for compiling .cpp files from the source directory (../)
$(OBJDIR)%.o: $(SRCDIR)%.cpp $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<

# Rule for compiling .cpp files from the library directory (../../lib)
$(OBJDIR)%.o: $(LIBDIR)%.cpp $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<

# --- DIRECTORY CREATION RULE (THE FIX) ---
# This rule tells make how to create the target directory.
# It's triggered by the `| $(TARGETDIR)` prerequisite in the main target rule.
# It uses Windows-specific commands.
$(TARGETDIR):
	@if not exist "$(subst /,\, $@)" mkdir "$(subst /,\, $@)"

# ==============================================================================
# --- CLEANUP ---
# ==============================================================================

# Rule to remove generated files
clean:
	-@del $(subst /,\,$(TARGET)) >nul 2>&1
	-@del *.o >nul 2>&1

# Phony targets are not files
.PHONY: all clean
