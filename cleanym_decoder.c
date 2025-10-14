#define PROGNAME "CLEANYM Decoder v1.0, by Mike Saarna. 2025."

#include <stdio.h>
#include <stdlib.h>
#include <sndfile.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <limits.h>

#include <unistd.h>

long TARGETRATE = 14000;

long samplecount;
float *samplebuffer = NULL;
FILE *infile;
SNDFILE *outfile;

int loadadpcm(char *filename, int preserveRate);
void usage(char *programname);
float getbitrate (long noisescore, long samplecount);

uint8_t encodesample_ym4(int16_t RawSample);
int16_t decodesample_ym4(uint8_t AdpcmSample);


// we keep these as globals, to easily save and restore for A:B testing
int16_t decoder_step_size = 127;
int32_t decoder_history = 0;

// step table used for 4-bit Decoder
const int step_table_4bit[8] = { 230, 230, 230, 230, 307, 409, 512, 614 }; 

int16_t save_decoder_step_size;
int32_t save_decoder_history;

int main(int argc, char **argv)
{
    int c;
    long int t;
    extern char *optarg;
    extern int optind, optopt;
    int errflg = 0;

    long insize;

    char *infilename, *outfilename;

    infilename = NULL;
    outfilename = NULL;

    while ((c = getopt(argc, argv, ":i:o:f:r:Rs")) != -1)
    {
	switch (c)
	{
	case 'i':
		infilename = optarg;
		break;
	case 'o':
		outfilename = optarg;
		break;
	case 'r':
		TARGETRATE = atoi(optarg);
		break;
	case ':': 	// option that should have operand, without operand
		fprintf(stderr, "*** ERR: option -%c requires an operand\n", optopt);
		errflg++;
		break;
	case '?':
		fprintf(stderr, "*** ERR: unrecognised option \"-%c\"\n", optopt);
		errflg++;
		break;
	}
    }
    if (errflg)
    {
	usage(argv[0]);
	return (2);
    }

    if (infilename == NULL)
    {
	fprintf(stderr, "-i wav file parameter required\n");
	usage(argv[0]);
	return 1;
    }

    if (outfilename == NULL)
    {
	fprintf(stderr, "Error: -o output file parameter required\n\n");
	usage(argv[0]);
	return 1;

    }

    infile=fopen(infilename,"rb");
    if (infile == NULL)
    {
	fprintf(stderr, "*** ERR: couldn't open \'%s\' for reading.\n",infilename);
	return 1;
    }

    SF_INFO sfinfo = {
        .samplerate = TARGETRATE,
        .channels = 1,          // Mono
        .format = SF_FORMAT_WAV | SF_FORMAT_FLOAT  // WAV, 32-bit float
    };


    outfile = sf_open(outfilename, SFM_WRITE, &sfinfo);
    if (outfile == NULL)
    {
        fprintf(stderr, "*** ERR: %s\n", sf_strerror(NULL));
	return 1;
    }

    fseek(infile,0L,SEEK_END);
    insize=ftell(infile);
    fseek(infile,0L,SEEK_SET);

    uint8_t *AdpcmBytes = calloc(insize,sizeof(uint8_t));
    if (AdpcmBytes == NULL)
    {
	fprintf(stderr, "*** ERR: couldn't allocate memory for input buffer.\n");
	return 1;
    }
    if(fread(AdpcmBytes,sizeof(uint8_t),insize,infile)!=insize)
    {
	fprintf(stderr, "*** ERR: Couldn't read %ld bytes from '%s'\n",insize,infilename);
	return 1;
    }
    fclose(infile);

    
    float *SampleBytes = calloc(insize*2,sizeof(float));
    if (SampleBytes == NULL)
    {
	fprintf(stderr, "*** ERR: couldn't allocate memory for output buffer.\n");
	return 1;
    }
   
    for (t=0;t<insize;t++)
    {
        SampleBytes[(t*2)+0]=((float)decodesample_ym4(AdpcmBytes[t]&0x0F))/32767.0;
        SampleBytes[(t*2)+1]=((float)decodesample_ym4((AdpcmBytes[t]>>4)&0x0F))/32767.0;
    }

    sf_writef_float(outfile, SampleBytes, insize*2);  
    
    sf_close(outfile);

    return (0);
}

#define CLAMP(x, low, high)  (((x) > (high)) ? (high) : (((x) < (low)) ? (low) : (x)))

int16_t decodesample_ym4(uint8_t AdpcmSample)
{
	// The 4 bit YMZ ADPCM Decoder
	// https://github.com/superctr/adpcm/blob/master/ymz_codec.c

	int32_t sign,delta,nstep,diff;

	// adjust step size and history
	sign  = AdpcmSample & 8;
	delta = AdpcmSample & 7;
	diff  = ((1+(delta<<1)) * decoder_step_size) >> 3;
	nstep = (step_table_4bit[delta] * decoder_step_size) >> 8;
	diff  = CLAMP(diff, 0, 32767);
	if (sign > 0)
		decoder_history -= diff;
	else
		decoder_history += diff;

	decoder_step_size = CLAMP(nstep, 127, 24576);
	decoder_history   = CLAMP(decoder_history, -32768, 32767);

	return(decoder_history);
}

void usage(char *programname)
{
    fprintf(stderr, "%s %s %s\n", PROGNAME, __DATE__, __TIME__);
    fprintf(stderr, "Usage: %s -i INPUTFILE -o OUTFILE [-r RATE]\n", programname);
    fprintf(stderr, "       INPUTFILE is a YMZ ADPCM file.\n");
    fprintf(stderr, "       -o specifies the output WAV file name.\n");
    fprintf(stderr, "       -r specifies the output bitrate in samples/second. (Default is %ld)\n",TARGETRATE);
    fprintf(stderr, "\n");
}

