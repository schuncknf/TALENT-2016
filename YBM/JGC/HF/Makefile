
#LIBS:=`root-config --libs`
#INCS:=`root-config --cflags`
ARM:= -Ofast

Esercizio: Esercizio.o oggetto.o
	g++ Esercizio.o oggetto.o -o Esercizio #${LIBS} ${INCS}

Esercizio.o: Esercizio.C
	g++ -c Esercizio.C -o Esercizio.o  #${LIBS} ${INCS}

oggetto.o: oggetto.C
	g++ -c oggetto.C -o oggetto.o  #${LIBS} ${INCS}

clean_ob:
	rm -rf *.o

clean: clean_ob
	rm -rf Esercizio

do: ;
	./Esercizio
