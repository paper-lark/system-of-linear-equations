.PHONY: all clean

SRC_DIR = src
BIN_DIR = build
DEBUG_FLAGS = -o ${BIN_DIR}/debug.out -g -fsanitize=bounds -fsanitize=address -fsanitize=undefined -lm
RELEASE_FLAGS = -o ${BIN_DIR}/release.out -O2 -lm

all: release debug

release: ${SRC_DIR}/main.c ${SRC_DIR}/functions.c ${SRC_DIR}/iteration.c ${SRC_DIR}/gauss.c
	gcc ${RELEASE_FLAGS} ${SRC_DIR}/main.c ${SRC_DIR}/functions.c ${SRC_DIR}/iteration.c ${SRC_DIR}/modified-gauss.c ${SRC_DIR}/gauss.c

debug: ${SRC_DIR}/main.c ${SRC_DIR}/functions.c ${SRC_DIR}/iteration.c ${SRC_DIR}/gauss.c
	gcc ${DEBUG_FLAGS} ${SRC_DIR}/main.c ${SRC_DIR}/functions.c ${SRC_DIR}/iteration.c ${SRC_DIR}/modified-gauss.c ${SRC_DIR}/gauss.c

clean:
	rm -rf ${BIN_DIR}/*