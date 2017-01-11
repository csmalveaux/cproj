##########################################################################
# buildMaze.mk
# Maze Builder Executable MakeFile
# Chloe Malveaux
# Dec 18 2016
# Version 1.0
##########################################################################

CD              = cd
MAKE            = make
RM              = rm -f
CP              = cp
CC			    = gcc

APP = MazeBuilder.exe

OUTPUTDIR   = obj
OBJDIR      = $(OUTPUTDIR)

INTOFILES   = $(MAZE_FILES)

MAZE_FILES = $(addprefix $(OUTPUTDIR)/, buildMaze.o)
MAZE_PATH  = ./

LDFLAGS = -lncurses

#CFLAGS = -fstack-protector-all

#LIBS = 

vpath	%.c $(MAZE_PATH)
vpath	%.h $(MAZE_PATH)

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
	$(CC) -c -g $(CFLAGS) $< -o $@

# Cleaning
.PHONY : clean ftp
clean :
	$(RM) $(INTOFILES) $(APP)