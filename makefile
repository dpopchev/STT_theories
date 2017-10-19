# used sources
#	https://www.gnu.org/software/make/manual/html_node/index.html#SEC_Contents
#	https://www.cs.swarthmore.edu/~newhall/unixhelp/howto_makefiles.html
#	https://stackoverflow.com/questions/9178285/how-can-makefile-use-separate-directories-for-source-code-and-binaries

# define the C compiler
CC = gcc

# define any compile-time flags
CFLAGS = -Wall -lm

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

# The following part of the makefile is combination of the generic parts of the sources listed above

.PHONY: all clean $(BIN_NAME) build run

all: $(BIN_NAME)
	@echo Simple compiler for simple needs

build:
	@echo We build objects
	$(CC) -c $(SRCS) -o $(OBJS) -I $(INCLUDES) $(CFLAGS)

$(BIN_NAME): build
	@echo We create binaries
	$(CC) $(OBJS) -o $(BIN_DIR)/$@

clean:
	rm -f $(BIN_DIR)/$(BIN_NAME) $(BUILD_DIR)/*.o

run:
	./$(BIN_DIR)/$(BIN_NAME)
