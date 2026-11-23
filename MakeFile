cc := g++
SRCS := main.cpp lageometry.cpp
OBJS := main.o lageometry.o
OUTPUTFILR := test
FLAGS := -g -O0 -lpthread

all: ${OBJS}
	${cc} $(FLAGS) main.o lageometry.o -o ${OUTPUTFILR}
	@echo ${OUTPUTFILR}
	@echo ${OBJS}
	@echo ${SRCS}

%.o: %.cpp
	@echo compling $<
	${cc} -c $(FLAGS) $< -o $@

install:
	make all

uninstall:
	make clean

clean:
	rm -f *.o
	rm -f ${OUTPUTFILR}
