CFLAGS=		-g -Wall -O2 -ftree-vectorize
DFLAGS=		-D_FILE_OFFSET_BITS=64 -D_USE_KNETFILE

indelsnpgen: indelsnpgen.o gencore.o indelsnpgenhelper.o seqreader.o libhts
	gcc $(CFLAGS) $(DFLAGS) indelsnpgen.o gencore.o indelsnpgenhelper.o seqreader.o \
	libhts.a /usr/lib/libgsl.a -lz -pthread -lstdc++ \
	-o indelsnpgen

indelsnpgen.o: indelsnpgen.cpp
	gcc $(CFLAGS) $(DFLAGS) -c indelsnpgen.cpp

gencore.o: gencore.cpp
	gcc $(CFLAGS) $(DFLAGS) -c gencore.cpp

indelsnpgenhelper.o: indelsnpgenhelper.cpp
	gcc $(CFLAGS) $(DFLAGS) -c indelsnpgenhelper.cpp

seqreader.o: seqreader.cpp
	gcc $(CFLAGS) $(DFLAGS) -c seqreader.cpp

libhts:
	make -C ~/htslib-0.2.0-rc10 libhts.a
	cp ~/htslib-0.2.0-rc10/libhts.a .

clean:
	rm -f indelsnpgen
	rm -f indelsnpgen.o
	rm -f gencore.o
	rm -f indelsnpgenhelper.o
	rm -f seqreader.o

