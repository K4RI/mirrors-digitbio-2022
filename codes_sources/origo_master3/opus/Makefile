VPATH = ../source
SRC = find_itemsets.cpp fisher.cpp globals.cpp itemset.cpp load_data.cpp opus_miner.cpp print_itemsets.cpp utils.cpp tidset.cpp itemQClass.cpp filter_itemsets.cpp find_closure.cpp
HDR = find_itemsets.h fisher.h globals.h itemset.h load_data.h opus_miner.h print_itemsets.h tidset.h utils.h itemQClass.h filter_itemsets.h find_closure.h
OBJ = ${SRC:.cpp=.o}
CFLAGS=-O3 -DNDEBUG
#CFLAGS=-pg -g -DNDEBUG

%.o : %.cpp $(HDR)
	gcc $(CFLAGS) -c $<

opus_miner: $(OBJ)
	g++ $(CFLAGS) -o $@ -lm $(OBJ)

opus_miner16: $(OBJ)
	g++ $(CFLAGS) -DSIXTEENBIT -o $@ -lm $(SRC)

opus_miner.tgz: $(SRC) $(HDR) Makefile README.txt
	tar -czf $@ $^
