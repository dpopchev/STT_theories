# used sources
#	https://www.gnu.org/software/make/manual/html_node/index.html#SEC_Contents
#	https://www.cs.swarthmore.edu/~newhall/unixhelp/howto_makefiles.html
#	https://stackoverflow.com/questions/9178285/how-can-makefile-use-separate-directories-for-source-code-and-binaries
#	https://stackoverflow.com/questions/7004702/how-can-i-create-a-makefile-for-c-projects-with-src-obj-and-bin-subdirectories
#	https://stackoverflow.com/questions/40621451/makefile-automatically-compile-all-c-files-keeping-o-files-in-separate-folde

# define the C compiler
CC = gcc

# define linker
LINKER = gcc

# default text editor
TEXT_EDITOR = mousepad

# default debugger / debug viewer
DEBUG_VIEW = nemiver

# define any compile-time flags
# used source for here
#	https://www.gnu.org/software/gsl/manual/html_node/GCC-warning-options-for-numerical-programs.html
# 	https://stackoverflow.com/questions/42586080/gcc-linking-object-files-with-warning-optimization-flags
CFLAGS_WARNS = -std=c11 -pedantic -Wall -W -Wmissing-prototypes -Wstrict-prototypes \
			   -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align \
			   -Wwrite-strings -Wnested-externs -Wextra -Wno-unused
CFLAGS_OPTIM = -fshort-enums -fno-common -Dinline= -march=native -O2
CFLAGS_LINKING = -lm -no-pie
CFLAGS_DEBUG = -g -std=c11

# define any directories containing header files other than /usr/include
INCLUDES_DIR = ./headers
INCLUDES = $(wildcard $(INCLUDES_DIR)/*.h)

# define the C source dirs and files
SRC_DIR = ./src
SRCS = $(wildcard $(SRC_DIR)/*.c)

# define place to put the object files and their names
# define the objects by replacing the .c of all words in the macro SRCS with the .o suffix
OBJ_DIR = ./obj
#OBJS = $(notdir $(SRCS:.c=.o))
OBJS = $(SRCS:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)

# define the executable file and where to be placed
BIN_DIR = ./bin
BIN_NAME = simulate
BIN_VALGRIND = simulate_valgrind
VALGRIND_TMP_FILE = qvasd
BIN_GDB = simulate_gdb

.PHONY: build
build: $(OBJS)
	@echo Objects have been created successfully

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c | $(OBJ_DIR)
	$(CC) $(CFLAGS_WARNS) $(CFLAGS_OPTIM) -I$(INCLUDES_DIR) -c $< -o $@

$(OBJ_DIR):
	@echo $@ folder missing, creating it
	mkdir -p $@

.PHONY: compile
compile: $(BIN_NAME)
	@echo Executable $(BIN_NAME) has been created successfully

$(BIN_NAME): build | $(BIN_DIR)
	@echo Creating executable $(BIN_NAME)
	$(LINKER) $(OBJS) -o $(BIN_DIR)/$@ $(CFLAGS_LINKING)

$(BIN_DIR):
	@echo $@ folder missing, creating it
	mkdir -p $@

.PHONY: run
run: | $(BIN_NAME)
	@echo Executing $(BIN_NAME)
	./src/python/live_plot_Results.py &
	./src/python/live_plot_solver_tries.py &
	./$(BIN_DIR)/$(BIN_NAME)

.PHONY: valgrind
valgrind: | $(BIN_DIR)
	@echo Will check for memory integrity using valgrind
	$(CC) -o $(BIN_DIR)/$(BIN_VALGRIND) $(SRCS) -I$(INCLUDES_DIR) \
		$(CFLAGS_DEBUG) $(CFLAGS_LINKING)
	valgrind -v --log-file="$(BIN_DIR)/$(VALGRIND_TMP_FILE)" \
		--read-var-info=yes --track-origins=yes --leak-check=full \
		--show-leak-kinds=all --tool=memcheck \
		./$(BIN_DIR)/$(BIN_VALGRIND)
	$(TEXT_EDITOR) $(BIN_DIR)/$(VALGRIND_TMP_FILE)

.PHONY: gdb
gdb: | $(BIN_DIR)
	@echo Will create gdb compatible executable
	$(CC) -o $(BIN_DIR)/$(BIN_GDB) $(SRCS) -I$(INCLUDES_DIR) \
		$(CFLAGS_DEBUG) $(CFLAGS_LINKING)
	$(DEBUG_VIEW) $(BIN_DIR)/$(BIN_GDB)

.PHONY: clean
clean:
	@echo Cleaning content of following directories: $(BIN_DIR), $(OBJ_DIR)
	rm -f \
		$(BIN_DIR)/$(BIN_NAME) \
		$(BIN_DIR)/$(BIN_VALGRIND) \
		$(BIN_DIR)/$(BIN_GDB) \
		$(BIN_DIR)/qvasd \
		$(OBJ_DIR)/*.o
