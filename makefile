# used sources
#	https://www.gnu.org/software/make/manual/html_node/index.html#SEC_Contents
#	https://www.cs.swarthmore.edu/~newhall/unixhelp/howto_makefiles.html
#	https://stackoverflow.com/questions/9178285/how-can-makefile-use-separate-directories-for-source-code-and-binaries

# define the C compiler
CC = gcc

# define any compile-time flags
# used source for here
#	https://www.gnu.org/software/gsl/manual/html_node/GCC-warning-options-for-numerical-programs.html
CFLAGS_WARNS = -std=c11 -pedantic -Wall -W -Wmissing-prototypes -Wstrict-prototypes \
			   -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align \
			   -Wwrite-strings -Wnested-externs -Wextra -Wno-unused
CFLAGS_OPTIM = -fshort-enums -fno-common -Dinline=
CFLAGS_LINKING = -lm
CFLAGS_DEBUG = -g -std=c11

# define any directories containing header files other than /usr/include
INCLUDES = ./headers

# define the C source dirs and files
SRC_DIR = ./src
SRCS = $(wildcard $(SRC_DIR)/*.c)

# define place to put the object files and their names
# define the objects by replacing the .c of all words in the macro SRCS with the .o suffix
BUILD_DIR = ./build
OBJS = $(BUILD_DIR)/$(notdir $(SRCS:.c=.o))

# define the executable file and where to be placed
BIN_DIR = ./bin
BIN_NAME = simulate
BIN_VALGRIND = simulate_valgrind
BIN_GDB = simulate_gdb

# The following part of the makefile is combination of the generic parts of the sources listed above
.PHONY: all clean $(BIN_NAME) build run valgrind_check gdb_tracking

all: $(BIN_NAME)
	@echo Simple compiler for simple needs

build:
	@echo Create build dir if necessary and create object files
	mkdir -p $(BUILD_DIR)
	$(CC) -c $(SRCS) -o $(OBJS) -I $(INCLUDES) $(CFLAGS_WARNS) $(CFLAGS_OPTIM)

$(BIN_NAME): build
	@echo Create bin dir if necessary and create executables in it
	mkdir -p $(BIN_DIR)
	$(CC) $(OBJS) -o $(BIN_DIR)/$@ $(CFLAGS_LINKING)

valgrind_check:
	@echo We will check the memory integrity of the code
	$(CC) -o $(BIN_DIR)/$(BIN_VALGRIND) $(SRCS) -I $(INCLUDES) $(CFLAGS_DEBUG) $(CFLAGS_LINKING)
	valgrind -v --log-file="$(BIN_DIR)/qvasd" \
		--read-var-info=yes --track-origins=yes --leak-check=full --show-leak-kinds=all \
		--tool=memcheck ./$(BIN_DIR)/$(BIN_VALGRIND)
	gvim $(BIN_DIR)/qvasd

gdb_tracking:
	@echo We will try to track the code with gdb
	$(CC) -o $(BIN_DIR)/$(BIN_GDB) $(SRCS) -I $(INCLUDES) $(CFLAGS_DEBUG) $(CFLAGS_LINKING)
	nemiver $(BIN_DIR)/$(BIN_GDB)

clean:
	rm -f \
		$(BIN_DIR)/$(BIN_NAME) \
		$(BIN_DIR)/$(BIN_VALGRIND) \
		$(BIN_DIR)/$(BIN_GDB) \
		$(BIN_DIR)/qvasd \
		$(BUILD_DIR)/*.o

run:
	./$(BIN_DIR)/$(BIN_NAME)
