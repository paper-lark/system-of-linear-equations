.PHONY: all clear

SRC_DIR = src
BIN_DIR = bin
DEBUG_FLAGS = -o ${BIN_DIR}/debug -g -fsanitize=bounds -fsanitize=address -fsanitize=undefined
RELEASE_FLAGS = -o ${BIN_DIR}/release -O2

all: release debug

release: ${SRC_DIR}/main.c ${SRC_DIR}/functions.c ${SRC_DIR}/iteration.c ${SRC_DIR}/gauss.c
	gcc ${RELEASE_FLAGS} ${SRC_DIR}/main.c ${SRC_DIR}/functions.c ${SRC_DIR}/iteration.c ${SRC_DIR}/modified-gauss.c ${SRC_DIR}/gauss.c

debug: ${SRC_DIR}/main.c ${SRC_DIR}/functions.c ${SRC_DIR}/iteration.c ${SRC_DIR}/gauss.c
	gcc ${DEBUG_FLAGS} ${SRC_DIR}/main.c ${SRC_DIR}/functions.c ${SRC_DIR}/iteration.c ${SRC_DIR}/modified-gauss.c ${SRC_DIR}/gauss.c

clear:
	rm *.o *.out