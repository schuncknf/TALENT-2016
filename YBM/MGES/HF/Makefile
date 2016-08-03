
#LIBS:=`root-config --libs`
#INCS:=`root-config --cflags`
ARM:= -O2

Esercizio: Esercizio.o oggetto.o #star.o
	g++ Esercizio.o oggetto.o -o Esercizio #-lncurses #${LIBS} ${INCS} 

Esercizio.o: Esercizio.C
	g++ -c Esercizio.C -o Esercizio.o -lncurses #${LIBS} ${INCS}

oggetto.o: oggetto.C
	g++ -c oggetto.C -o oggetto.o -lncurses #${LIBS} ${INCS}

#star.o: star.C
#	g++ -c star.C -o star.o -lncurses #${LIBS} ${INCS}

clean_ob:
	rm -rf *.o

clean: clean_ob
	rm -rf Esercizio

do: ;
	./Esercizio
