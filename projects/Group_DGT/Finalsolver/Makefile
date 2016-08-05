FF = gfortran
EXEC = mainfinal

CFLAGS = 
LDFLAGS =


objects = variables.o potentialsfinal.o calcwavesfinal.o mainfinal.o

.SUFFIXES: .f90 .o


$(EXEC): $(objects)
	$(FF) -o $(EXEC) $(objects) $(LDFLAGS)

.f90.o:
	$(FF) $(CFLAGS) -o $(@) $<;


modules:
	$(FF) $(CFLAGS) $(mods) 

%.o: %.f90
	$(FF) -o $@ -c $^ $(CFLAGS) 


clean:
	rm -rf *.o

mrproper: clean
	rm -rf $(EXEC)
