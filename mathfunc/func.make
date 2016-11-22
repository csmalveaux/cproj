###################################################################
#
# Copyright (c) 2001
# Boeing-SVS, Inc.
# 4411 The 25 Way NE, Suite 350
# Albuquerque, NM 87109
# All Rights Reserved
#
# Description:
#
# -I"C:\cygwin64\usr\include" -L"C:\cygwin64\usr\lib" -lgslcblas -lgsl -lm -o test

###################################################################


CD              = cd
MAKE            = make
RM              = rm -f
CP              = cp
CC			= gcc

CYGPATH		= C:/cygwin
CDIR		= C:/Data/Workspaces/ccode

APP			= mathfunc.o

INTOFILES		= $(MATHFUNC_FILES)

MATHFUNC_FILES = $(addprefix $(OBJDIR)/, test_functions.o sort.o search.o array.o binarytree.o matrix.o matrixDecomp.o)
MATHFUNC_PATH  = ./

LDFLAGS = -lm -I$(MATHFUNC_PATH) 

vpath	%.c $(MATHFUNC_PATH)
vpath	%.h $(MATHFUNC_PATH)

###################################################################

# Make rules
###################################################################

# All Definition
all:	$(APP)
	
# Linking
$(APP) :	$(INTOFILES)
	$(CC) -o $(APP) $(INTOFILES) $(LDFLAGS) $(LIBS)
	

# Handle Directory
$(INTOFILES) : | $(OBJDIR) $(OUTPUTDIR)
$(OBJDIR) :
	mkdir $(OBJDIR)
$(OUTPUTDIR) :
	mkdir $(OUTPUTDIR)
	
# Compiling
$(OBJDIR)/%.o : %.c
	$(CC) -c $(LDFLAGS) $< -o $@

# Cleaning
.PHONY : clean ftp
clean :
	$(RM) $(INTOFILES) $(APP)
