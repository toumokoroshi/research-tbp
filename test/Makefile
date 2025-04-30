CC = g++
CFLAGS = -Wall -Wextra -Wshadow -fno-stack-protector -static-libgcc -fopenmp -static-libstdc++ -std=c++20 

TARGETDIR = ./
TARGET = test.exe
SRCDIR = ./
LIBDIR = ../lib
OBJDIR = ./

BOOSTDIR = C:/BOOST/boost_1_87_0
INCLUDES = -I$(LIBDIR)
INCLUDES += -I$(BOOSTDIR)

 
SRCS = $(SRCDIR)test.cpp  
OBJS = $(OBJDIR)test.o 


HEADERS = $(LIBDIR)/crtbp.hpp $(LIBDIR)/vector3d.hpp $(LIBDIR)/utils.hpp  $(LIBDIR)/chaosmap_generator.hpp  $(LIBDIR)/rtbp.hpp

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