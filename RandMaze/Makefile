
###################################################################
#
###################################################################


CD              = cd
MAKE            = make
RM              = rm -f
CP              = cp
CC			    = gcc

APP			= test_maze.o

OUTPUTDIR   = obj
OBJDIR      = $(OUTPUTDIR)

INTOFILES		= $(MAZE_FILES) $(MATHFUNC_FILES)

MATHFUNC_FILES = $(addprefix $(OBJDIR)/, complex.o dsp.o matrix.o array.o singValD.o matrixDecomp.o sort.o search.o binarytree.o)
MATHFUNC_PATH  = $../../mathfunc

MAZE_FILES = $(addprefix $(OBJDIR)/, maze.o mazeDisp.o mazeGen.o mouse.o mazeSolv.o main.o)
MAZE_PATH  = ./

LDFLAGS = -I$(MATHFUNC_PATH)

CFLAGS = -fstack-protector-all

vpath	%.c $(MAZE_PATH):$(MATHFUNC_PATH)
vpath	%.h $(MAZE_PATH):$(MATHFUNC_PATH)

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
	$(CC) -c -g $(LDFLAGS) $(CFLAGS) $< -o $@

# Cleaning
.PHONY : clean ftp
clean :
	$(RM) $(INTOFILES) $(APP)