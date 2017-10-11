objects = main.o

my_program : $(objects)
	gcc -o my_program $(objects)

main.o :

.PHONY : clean
clean :
	rm my_program $(objects)

