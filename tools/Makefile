convert : plfs_xtc.o convert.o  \
	gcc -o convert convert.o plfs_xtc.o -lm

plfs_xtc.o : plfs_xtc.c plfs_xtc.h  \
	gcc -c plfs_xtc.c

convert.o : convert.c  \
	gcc -c convert.c

.PHONE : clean
clean : 
	-rm convert *.o 