# used sources
#	https://www.gnu.org/software/make/manual/html_node/index.html#SEC_Contents
#	https://www.cs.swarthmore.edu/~newhall/unixhelp/howto_makefiles.html

# define the C compiler
CC = gcc

# define any compile-time flags
CFLAGS = -Wall

# define any directories containing header files other than /usr/include
#INCLUDES =

# define the C source files
SRCS = main.c

# define the objects by replacing the .c of all words in the macro SRCS with the .o suffix
OBJS = $(SRCS:.c=.o)

# define the executable file
MAIN = simulate

#
# The following part of the makefile is generic; it can be used to
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#

.PHONY: depend clean

all: $(MAIN)
	@echo  Simple compiler for simple needs

$(MAIN): $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file)
# (see the gnu make manual section about automatic variables)
.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) *.o *~ $(MAIN)

depend: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it
