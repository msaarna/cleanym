# The cleanym generic makefile

CC     = gcc
LIBS   = -lsndfile -lsamplerate -lm
DEBUG  = -Wall -g
 
all: cleanym_encoder cleanym_decoder
 
cleanym_encoder:	cleanym_encoder.c
	$(CC) cleanym_encoder.c $(DEBUG) $(LIBS) -o cleanym_encoder

cleanym_decoder:	cleanym_decoder.c
	$(CC) cleanym_decoder.c $(DEBUG) $(LIBS) -o cleanym_decoder

clean:
	rm -f *.o cleanym_encoder cleanym_decoder

backup:
	rm -f samples/2_* samples/3_*
	cd .. ; tar -cvzf cleanym_encoder.tgz cleanym_encoder
	cp ../cleanym_encoder.tgz /volume1/data/public/_Mike/backups/
