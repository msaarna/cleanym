#define PROGNAME "CLEANYM Encoder v1.0, by Mike Saarna. 2025."

#include <stdio.h>
#include <stdlib.h>
#include <sndfile.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <float.h>

#include <samplerate.h>

#define SYNC_QUALITY SRC_SINC_MEDIUM_QUALITY

#define LOOKAHEAD_DEPTH 4
// ^- The number of samples to consider for lookahead optimisation.
// 4 samples seems to be the sweet spot for quality and speed, for 2 reasons: 
//   1. Adding 5+ samples adds weight to lower-frequency errors, but higher-
//      frequency ADPCM errors are always the most glaring.
//   2. Due to look-ahead changes we may subsequently perform, the further 
//      into the future we look, the less accurate the data is.

// Weights for the enhanced psychoacoustic error metric.
#define ERROR_WEIGHT_ABSOLUTE   1.00f
#define ERROR_WEIGHT_VOLATILITY 2.00f 

long TARGETRATE = 32160;

#include <unistd.h>

#define FRAMESIZE 64
long SAMPLERATE;

enum outputformats
{ RAWFMT = 0, ASMFMT = 1, CFMT = 2 };
int outputformat = RAWFMT;

int64_t samplecount;
float *samplebuffer = NULL;
FILE *outfile;

int loadwave (char *filename, int preserveRate);
void usage (char *programname);
float getbitrate (long noisescore, int64_t samplecount);
float normalizeSample (float *normBuffer, int64_t normBufferSize);
void apply_hipass_filter_to_buffer (float *buffer, int num_samples, float cutoff_freq_hz, float sample_rate_hz);
float returnMedianValue (float a, float b, float c, float d, float e);
int compare_floats (const void *a, const void *b);

uint8_t encodesample_ym4 (int16_t RawSample);
int16_t decodesample_ym4 (uint8_t AdpcmSample);
uint8_t encodesample_ym4_lookahead (int64_t current_sample_index, int16_t CurrentRawSample, long lookbehind_error);

void save_all (void);
void save_encoder_state (void);
void save_decoder_state (void);
void restore_all (void);
void restore_encoder_state (void);
void restore_decoder_state (void);

// used for stats
long noisescore = 0;
long framenoise = 0;
long worstframenoise = 0;
long fixcount = 0;

// used for normalisation
int8_t volumeDivisor = 1;

// Max slew rate default. Smaller values = more aggressive limiting.
float max_slew_rate_factor = 0.75f;

// step table used for 4-bit Encoder and Decoder
const int step_table_4bit[8] = { 230, 230, 230, 230, 307, 409, 512, 614 };

// we keep these as globals, to easily save and restore for A:B testing
int16_t encoder_step_size = 127;
int32_t encoder_history = 0;
int16_t decoder_step_size = 127;
int32_t decoder_history = 0;

int16_t save_encoder_step_size;
int32_t save_encoder_history;
int16_t save_decoder_step_size;
int32_t save_decoder_history;

// Define M_PI if not available (e.g., on MSVC before C++20 or in strict C modes)
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Hi-pass filter state 
static float hp_x_prev = 0.0f;	// Previous input sample x[n-1]
static float hp_y_prev = 0.0f;	// Previous output sample y[n-1]
static float hp_alpha = 0.0f;	// Filter coefficient, calculated once

int main (int argc, char **argv)
{
    int c;
    int64_t t;
    extern char *optarg;
    extern int optind, optopt;
    int errflg = 0;

    int preserveRate = 0;

    char *infilename, *outfilename;

    infilename = NULL;
    outfilename = NULL;

    while ((c = getopt (argc, argv, ":i:o:f:r:R")) != -1)	// getopt string unchanged
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
	    TARGETRATE = atoi (optarg);
	    break;
	case 'R':
	    preserveRate = 1;
	    break;
	case ':':
	    fprintf (stderr, "*** ERR: option -%c requires an operand\n", optopt);
	    errflg++;
	    break;
	case '?':
	    fprintf (stderr, "*** ERR: unrecognised option \"-%c\"\n", optopt);
	    errflg++;
	    break;
	}
    }
    if (errflg)
    {
	usage (argv[0]);
	return (2);
    }

    if (infilename == NULL)
    {
	fprintf (stderr, "-i wav file parameter required\n\n");
	usage (argv[0]);
	return 1;
    }

    if (outfilename == NULL)
    {
	fprintf (stderr, "-o output file parameter required\n\n");
	usage (argv[0]);
	return 1;

    }

    printf ("Loading file %s\n", infilename);
    if (loadwave (infilename, preserveRate) != 0)
    {
	fprintf (stderr, "error: loading wav file failed\n");
	// Ensure samplebuffer is freed if loadwave fails after allocation
	if (samplebuffer != NULL)
	{
	    free (samplebuffer);
	    samplebuffer = NULL;
	}
	return 1;
    }

    // At this point the input WAV file has been loaded and converted 
    // into floats, ranging from -32512 to 32512

    outfile = fopen (outfilename, "wb");
    if (outfile == NULL)
    {
	fprintf (stderr, "error: creation of output file \"%s\" failed\n", outfilename);
	if (samplebuffer != NULL)
	    free (samplebuffer);
	return 1;
    }
    strcat (outfilename,".raw");
    FILE *outfile2 = fopen(outfilename,"wb");

    uint8_t AdpcmSample = 0;
    uint32_t AdpcmSamplePacked = 0;
    int16_t SourcePcmSample, PcmSample;
    long last_actual_error = 0; // For look-behind error tracking.

    float *frame_distortions = calloc ((samplecount / FRAMESIZE) + 1, sizeof (float));
    if (frame_distortions == NULL)
    {
	fprintf (stderr, "Could not allocate memory for tracking distortions.\n");
	fclose (outfile);
	if (samplebuffer != NULL)
	    free (samplebuffer);
	return 1;
    }

    for (t = 0; t < samplecount; t++)
    {

	if ((t % FRAMESIZE) == 0 && t > 0)
	{
	    // end of a frame - record info for noise stats...
	    if (framenoise > worstframenoise)
		worstframenoise = framenoise;

	    // Ensure index is within bounds for frame_distortions
	    if ((t / FRAMESIZE) < (samplecount / FRAMESIZE) + 1)
	    {			// Check should be against allocated size
		frame_distortions[t / FRAMESIZE] = (float) framenoise / FRAMESIZE;
	    }

	    framenoise = 0;
	}

	SourcePcmSample = (int16_t) samplebuffer[t];	// Fetch this wav sample (already scaled to int16 range as float).

	// Use the lookahead encoder if there's at least one next sample.
	// The lookahead function itself handles cases where fewer than (LOOKAHEAD_DEPTH-1) future samples are available.
	if (t < samplecount - 1)
	{	// If there's at least one next sample for lookahead
        // Pass the last actual error to the lookahead function.
	    AdpcmSample = encodesample_ym4_lookahead (t, SourcePcmSample, last_actual_error);
	}
	else
	{	// For the very last sample, no lookahead is possible with this scheme; use normal encoder.
	    AdpcmSample = encodesample_ym4 (SourcePcmSample);
	}


	// Run the decoder, just for the upcoming noise stats.
	// The encodesample_ym4_lookahead function ensures that global decoder_step_size/history
	// are restored to their state *before* this current sample's AdpcmSample was determined.
	// So, this decodesample_ym4 call uses the correct initial decoder state and advances it,
	// keeping it synchronized with the encoder's state progression.
	PcmSample = (int16_t) decodesample_ym4 (AdpcmSample);

    // Calculate and store the actual error for the next iteration's look-behind.
    last_actual_error = (long)PcmSample - (long)SourcePcmSample;

	noisescore += abs (PcmSample - SourcePcmSample);
	framenoise += abs (PcmSample - SourcePcmSample);

	fputc(((PcmSample - SourcePcmSample)/256),outfile2);

	if (!(t & 1))		// first nibble
	{
	    AdpcmSamplePacked = (AdpcmSample & 0x0F);
	    // Make sure to handle the last sample if samplecount is odd
	    if (t == samplecount - 1)
	    {
		// If it's the last sample and it's an odd one
	        fputc ((char) AdpcmSamplePacked & 0xFF, outfile);
	    }
	    else
	    {
		continue;	// run the loop for the second nibble
	    }
	}
	else			// second nibble
	{
	    AdpcmSamplePacked |= (AdpcmSample & 0x0F) << 4;
	}

	if (t & 1)
	{
	    fputc ((char) AdpcmSamplePacked & 0xFF, outfile);
	}
	else if (t == samplecount - 1)
	{
	    // Handle the very last nibble if samplecount is odd
	    // This means only AdpcmSamplePacked = (AdpcmSample & 0x0F) was set.
	    // The byte should contain this nibble, other nibble can be 0.
	    fputc ((char) AdpcmSamplePacked & 0xFF, outfile);	// high nibble will be 0
	}

    }				// End of for loop over samples

    fclose (outfile);
    fclose (outfile2);

    printf ("Lookahead optimisation adjusted %ld 4-bit samples. (%.2f%%)\n", fixcount, ((float)fixcount/(float)samplecount)*100.0f);
    printf ("Wrote %lld 4-bit samples, packed into %lld bytes.\n", (long long) samplecount, (long long) (samplecount + 1) / 2);

    if (samplecount / FRAMESIZE > 0)
    {
	qsort (frame_distortions, samplecount / FRAMESIZE, sizeof (float), compare_floats);
    }

    float avgbits, worstbits, medianbits;

    avgbits = getbitrate (noisescore, samplecount);
    worstbits = getbitrate (worstframenoise, (FRAMESIZE > 0 ? FRAMESIZE : 1));

    printf ("\n");
    printf ("ADPCM quality report...\n");

    if (samplecount / FRAMESIZE > 0)
    {
	medianbits = getbitrate ((long) frame_distortions[(samplecount / FRAMESIZE) / 2], 1);	// Median distortion is per-sample average for that frame
	printf ("  avg    distortion: %*ld    avg bits/sample:    %2.2f\n", 6,
		(samplecount > 0 ? noisescore / samplecount : 0L), avgbits);
	printf ("  median distortion: %*ld    median bits/sample: %2.2f\n", 6,
		(long) frame_distortions[(samplecount / FRAMESIZE) / 2], medianbits);
	printf ("  worst  distortion: %*ld    worst bits/sample:  %2.2f\n", 6,
		(FRAMESIZE > 0 ? worstframenoise / FRAMESIZE : 0L), worstbits);
    }
    else
    {
	medianbits = 0.0f;
	printf ("  avg    distortion: %*ld    avg bits/sample:    %2.2f\n", 6,
		(samplecount > 0 ? noisescore / samplecount : 0L), avgbits);
	printf ("  (median distortion not available for very short samples)\n");
	printf ("  worst  distortion: %*ld    worst bits/sample:  %2.2f\n", 6,
		(FRAMESIZE > 0 && worstframenoise > 0 ? worstframenoise / FRAMESIZE : 0L), worstbits);
    }
    printf ("\n");
    printf ("  Note: ADPCM bits/sample aren't directly comparable to PCM bits/sample.\n");
    printf ("  ADPCM noise increases as volume changes, masking it's perception.\n");
    

    free (frame_distortions);
    if (samplebuffer != NULL)
    {
	free (samplebuffer);
	samplebuffer = NULL;
    }
    return (0);
}

// Function signature and implementation to use new scoring metric.
uint8_t encodesample_ym4_lookahead (int64_t current_sample_index, int16_t CurrentRawSample, long lookbehind_error)
{
    uint8_t best_adpcm_for_current_sample = 0;
    float min_total_score = FLT_MAX;

    // Snapshot of global states at function entry.
    int16_t entry_encoder_step_size = encoder_step_size;
    int32_t entry_encoder_history = encoder_history;
    int16_t entry_decoder_step_size = decoder_step_size;
    int32_t entry_decoder_history = decoder_history;

    // For fixcount statistics: determine what the normal encoder would have chosen.
    save_all();
    uint8_t normal_adpcm_for_current_sample = encodesample_ym4(CurrentRawSample);
    restore_all();

    uint8_t trial_adpcm_current;
    for (trial_adpcm_current = 0; trial_adpcm_current <= 15; ++trial_adpcm_current)
    {
        // Store errors for the entire path for later scoring.
        long path_errors[LOOKAHEAD_DEPTH] = {0};
        int path_length = 0;

        int16_t sim_path_encoder_step_size;
        int32_t sim_path_encoder_history;

        // --- Part 1: Process CurrentRawSample with trial_adpcm_current ---
        decoder_step_size = entry_encoder_step_size;
        decoder_history = entry_encoder_history;

        int16_t decoded_pcm_current = decodesample_ym4(trial_adpcm_current);
        path_errors[0] = (long)decoded_pcm_current - (long)CurrentRawSample;
        path_length++;

        sim_path_encoder_step_size = decoder_step_size;
        sim_path_encoder_history = decoder_history;

        // --- Part 2: Process subsequent (LOOKAHEAD_DEPTH - 1) samples ---
        int i;
        for (i = 1; i < LOOKAHEAD_DEPTH; ++i)
        {
            int64_t future_sample_index = current_sample_index + i;
            if (future_sample_index >= samplecount)
            {
                break;
            }
            int16_t FutureRawSample = (int16_t)samplebuffer[future_sample_index];

            encoder_step_size = sim_path_encoder_step_size;
            encoder_history = sim_path_encoder_history;
            uint8_t adpcm_for_future_sample = encodesample_ym4(FutureRawSample);

            decoder_step_size = sim_path_encoder_step_size;
            decoder_history = sim_path_encoder_history;
            int16_t decoded_pcm_future = decodesample_ym4(adpcm_for_future_sample);
            
            path_errors[i] = (long)decoded_pcm_future - (long)FutureRawSample;
            path_length++;

            sim_path_encoder_step_size = encoder_step_size;
            sim_path_encoder_history = encoder_history;
        }

        // --- Part 3: Calculate weighted psychoacoustic score ---
        float absolute_error_sum = 0.0f;
        float volatility_error_sum = 0.0f;
        long previous_error = lookbehind_error;

        for (i = 0; i < path_length; i++)
        {
            absolute_error_sum += fabsf((float)path_errors[i]);
            volatility_error_sum += fabsf((float)path_errors[i] - (float)previous_error);
            previous_error = path_errors[i];
        }
        
        float current_total_score = (ERROR_WEIGHT_ABSOLUTE * absolute_error_sum) + (ERROR_WEIGHT_VOLATILITY * volatility_error_sum);

        // --- Part 4: Compare total score and update best choice ---
        if (current_total_score < min_total_score)
        {
            min_total_score = current_total_score;
            best_adpcm_for_current_sample = trial_adpcm_current;
        }
    }

    // --- Finalize: Apply the best choice and set global states correctly ---
    if (normal_adpcm_for_current_sample != best_adpcm_for_current_sample)
        fixcount++;

    // Set the *actual global encoder state* for the chosen path.
    decoder_step_size = entry_encoder_step_size;
    decoder_history = entry_encoder_history;
    decodesample_ym4(best_adpcm_for_current_sample); 

    encoder_step_size = decoder_step_size;
    encoder_history = decoder_history;

    // Restore the *global decoder state* to what it was at function entry for the main loop's stats.
    decoder_step_size = entry_decoder_step_size;
    decoder_history = entry_decoder_history;

    return best_adpcm_for_current_sample;
}


float getbitrate (long noisescore, int64_t samplecount)
{
    // Ensure samplecount is not zero to prevent division by zero
    if (samplecount == 0)
    {
	return 0.0f;
    }

    long avgdistortion_long = noisescore / samplecount;

    // Handle case where average distortion is 0 (perfect encoding for this block)
    if (avgdistortion_long <= 0)
    {				// Can be 0 if noisescore is 0 or noisescore < samplecount
	return 16.0f;		// No bits lost to noise
    }

    float avgdistortion_float = (float) avgdistortion_long;

    // Calculate the number of bits needed to represent this average distortion.
    // This is log base 2 of the average distortion.
    // log2(x) = log(x) / log(2) using natural logarithm (log) or base-10 (log10)
    float bits_lost_to_noise = logf (avgdistortion_float) / logf (2.0f);

    // If average error is less than 1 LSB (in 16-bit context)
    // logf of value < 1 is negative. This means bits_lost_to_noise would be negative.
    // Consider such small errors as effectively zero bits lost.
    if (avgdistortion_float < 1.0f)
    {
	bits_lost_to_noise = 0.0f;
    }

    float effective_bitrate = 16.0f - bits_lost_to_noise;

    // Clamp effective_bitrate
    if (effective_bitrate < 0.0f)
    {
	effective_bitrate = 0.0f;
    }
    if (effective_bitrate > 16.0f) // highly unlikely
    {
	effective_bitrate = 16.0f;
    }

    return effective_bitrate;
}

float normalizeSample (float *normBuffer, int64_t normBufferSize)
{
    int64_t t;
    float maxvolume = 0.0f;
    float fVolumeDivisor;

    if (normBufferSize == 0)
	return 1.0f;		// Handle empty buffer case

    for (t = 0; t < normBufferSize; t++)
	if (fabsf (normBuffer[t]) > maxvolume)
	    maxvolume = fabsf (normBuffer[t]);

    if (maxvolume == 0.0f)	// Completely silent sample
    {
	fVolumeDivisor = 1.0f;	// No change needed
    }
    else if (maxvolume > 1.0f)	// Needs attenuation
    {
	// Attenuate to make peak 1.0f
	float attenuation_factor = 1.0f / maxvolume;
	for (t = 0; t < normBufferSize; t++)
	    normBuffer[t] = normBuffer[t] * attenuation_factor;
	fVolumeDivisor = 1.0f;
    }
    else			// Needs amplification (maxvolume is > 0 and <= 1.0)
    {
	// Amplify by a power of 2 to make peak close to 1.0f without exceeding it.
	fVolumeDivisor = exp2f (floorf (log2f (1.0f / maxvolume)));

	// Safety clamp for fVolumeDivisor, considering it's cast to int8_t later.
	// A very small maxvolume could lead to a huge fVolumeDivisor.
	if (fVolumeDivisor > 127.0f)
	    fVolumeDivisor = 127.0f;
	if (fVolumeDivisor < 1.0f)
	    fVolumeDivisor = 1.0f;	// Should not happen if maxvolume <= 1.0f and > 0

	for (t = 0; t < normBufferSize; t++)
	    normBuffer[t] = normBuffer[t] * fVolumeDivisor;
    }
    return fVolumeDivisor;
}

int compare_floats (const void *a, const void *b)
{
    float fa = *(const float *) a;
    float fb = *(const float *) b;
    if (fa < fb)
	return -1;
    if (fa > fb)
	return 1;
    return 0;
}

#define CLAMP(x, low, high)  (((x) > (high)) ? (high) : (((x) < (low)) ? (low) : (x)))

uint8_t encodesample_ym4 (int16_t RawSample)
{
    // The 4 bit YMZ ADPCM Encoder from here...
    // https://github.com/superctr/adpcm/blob/master/ymz_codec.c

    int32_t sign, step, delta, nstep, diff;
    unsigned int AdpcmSampleVal;

    step = RawSample - encoder_history;

    if (encoder_step_size == 0)
    {				// Prevent division by zero, though CLAMP should keep it >= 127
	AdpcmSampleVal = (abs (step) == 0) ? 0 : 7;
    }
    else
    {
	AdpcmSampleVal = (abs (step) << 16) / (encoder_step_size << 14);
    }
    AdpcmSampleVal = CLAMP (AdpcmSampleVal, 0, 7);

    if (step < 0)
	AdpcmSampleVal = AdpcmSampleVal | 8;

    // adjust step size and history
    sign = AdpcmSampleVal & 8;
    delta = AdpcmSampleVal & 7;

    diff = ((1 + (delta << 1)) * encoder_step_size) >> 3;
    nstep = (step_table_4bit[delta] * encoder_step_size) >> 8;
    diff = CLAMP (diff, 0, 32767);	// Max diff can be 46080, so this clamp is important.
    if (sign > 0)
	encoder_history -= diff;
    else
	encoder_history += diff;

    encoder_step_size = CLAMP (nstep, 127, 24576);
    encoder_history = CLAMP (encoder_history, -32768, 32767);

    return (AdpcmSampleVal & 0x0F);	// Ensure it's a nibble
}

int16_t decodesample_ym4 (uint8_t AdpcmSample)
{
    // The 4 bit YMZ ADPCM Decoder from here...
    // https://github.com/superctr/adpcm/blob/master/ymz_codec.c

    int32_t sign, delta, nstep, diff;

    sign = AdpcmSample & 8;
    delta = AdpcmSample & 7;

    diff = ((1 + (delta << 1)) * decoder_step_size) >> 3;
    nstep = (step_table_4bit[delta] * decoder_step_size) >> 8;
    diff = CLAMP (diff, 0, 32767);	// Consistent with encoder, max diff can be 46080.
    if (sign > 0)
	decoder_history -= diff;
    else
	decoder_history += diff;

    decoder_step_size = CLAMP (nstep, 127, 24576);
    decoder_history = CLAMP (decoder_history, -32768, 32767);

    return ((int16_t) decoder_history);
}

int loadwave (char *filename, int preserveRate)
{
    SF_INFO sndInfo;
    SNDFILE *sndFile;
    long numFrames;
    int channels;
    float *samplebuffertmp;
    int64_t t_load;

    // Open sound file
    sndFile = sf_open (filename, SFM_READ, &sndInfo);
    if (sndFile == NULL)
    {
	fprintf (stderr, "Error reading wav file '%s': %s\n", filename, sf_strerror (sndFile));
	return 1;
    }

    channels = sndInfo.channels;

    SAMPLERATE = sndInfo.samplerate;
    if (preserveRate)
	TARGETRATE = SAMPLERATE;

    samplecount = sndInfo.frames;

    samplebuffer = calloc (samplecount, sizeof (float));	// Allocate for actual number of frames for single channel float
    if (samplebuffer == NULL)
    {
	fprintf (stderr, "Could not allocate memory for file (samplebuffer)\n");
	sf_close (sndFile);
	return 1;
    }

    // Load the sample data into a array of floats...
    if (channels == 1)
	numFrames = sf_readf_float (sndFile, samplebuffer, samplecount);	// Use samplecount (int64_t)
    else			// channels>1
    {
	// we have at least two channels, but we only want one.
	// drop all channels except the first (left) one...
	long s = 0;
	samplebuffertmp = calloc (samplecount * channels, sizeof (float));	// Allocate for all channels' data
	if (samplebuffertmp == NULL)
	{
	    if (samplebuffer != NULL)
		free (samplebuffer);
	    samplebuffer = NULL;	// Mark as freed
	    fprintf (stderr, "Could not allocate memory for file (samplebuffertmp)\n");
	    sf_close (sndFile);
	    return 1;
	}

	numFrames = sf_readf_float (sndFile, samplebuffertmp, samplecount);

	for (t_load = 0; t_load < numFrames; t_load++)
	{
	    samplebuffer[t_load] = samplebuffertmp[s];
	    s = s + channels;
	}
	free (samplebuffertmp);
    }

    // Check correct number of samples loaded
    if (numFrames != samplecount)	// Compare with samplecount (int64_t)
    {
	fprintf (stderr, "Did not read enough frames for source. Expected %lld, got %ld\n", (long long) samplecount,
		 numFrames);
	sf_close (sndFile);
	if (samplebuffer != NULL)
	    free (samplebuffer);
	samplebuffer = NULL;
	return 1;
    }
    sf_close (sndFile);

    if (TARGETRATE == SAMPLERATE)
	printf ("The input samplerate is equal to output sample rate.\nNo scaling required\n");
    else
    {
	double ratio = (double) TARGETRATE / SAMPLERATE;
	int64_t input_frames = samplecount;
	// Calculate output_frames carefully to avoid overflow if samplecount is large.
	int64_t output_frames = (int64_t) round ((double) samplecount * ratio);

	float *samplebufferout = calloc (output_frames, sizeof (float));
	if (samplebufferout == NULL)
	{
	    fprintf (stderr, "Could not allocate memory for resampled buffer\n");
	    if (samplebuffer != NULL)
		free (samplebuffer);
	    samplebuffer = NULL;
	    return 2;		// Or appropriate error code
	}


	printf ("Scaling sample rate from %ld to %ld...", SAMPLERATE, TARGETRATE);
	fflush (stdout);

	SRC_DATA src_data = {
	    .data_in = samplebuffer,
	    .input_frames = input_frames,
	    .data_out = samplebufferout,
	    .output_frames = output_frames,
	    .src_ratio = ratio
	};

	int error = src_simple (&src_data, SYNC_QUALITY, 1);	// 1 channel
	if (error)
	{
	    fprintf (stderr, "\nResampling error: %s\n", src_strerror (error));
	    free (samplebufferout);
	    if (samplebuffer != NULL)
		free (samplebuffer);
	    samplebuffer = NULL;
	    return 2;
	}

	free (samplebuffer);

	samplebuffer = samplebufferout;
	samplecount = src_data.output_frames_gen;


	printf (" done. Output frames: %lld\n", (long long) samplecount);
    }

    // In case we encountered overly-hot samples...
    float norm_factor = normalizeSample (samplebuffer, samplecount);
    volumeDivisor = (int8_t) roundf (norm_factor);	// roundf and cast to int8_t
    if (volumeDivisor == 0 && norm_factor > 0)
	volumeDivisor = 1;	// Avoid divisor being 0

    fprintf (stderr, "Applied temporary normalization (factor: %f, stored divisor: %d).\n", norm_factor, volumeDivisor);


    // In case we introduced some issues with processing, re-normalize to [-1,1]
    // and get amplification factor if any.
    normalizeSample (samplebuffer, samplecount);	// This normalizes to [-1,1] if needed.

    // scale down by the factor determined by the *first* normalization, if it was an amplification.
    if (volumeDivisor > 1)
    {				// Only apply if original normalization decided to amplify
	for (t_load = 0; t_load < samplecount; t_load++)
	    samplebuffer[t_load] = samplebuffer[t_load] / volumeDivisor;
    }

    // Scale up to 16-bit signed integer range for ADPCM encoding
    for (t_load = 0; t_load < samplecount; t_load++)
	samplebuffer[t_load] = CLAMP (samplebuffer[t_load] * 32767.0f, -32768.0f, 32767.0f);

    return 0;
}

void hipass_filter_setup (float cutoff_freq_hz, float sample_rate_hz)
{
    if (cutoff_freq_hz <= 0.0f || sample_rate_hz <= 0.0f || cutoff_freq_hz >= sample_rate_hz / 2.0f)
    {
	hp_alpha = 1.0f;	// Pass-through or disable if params are invalid/extreme
	hp_x_prev = 0.0f;
	hp_y_prev = 0.0f;
	return;
    }
    double rc = 1.0 / (2.0 * M_PI * cutoff_freq_hz);
    double dt = 1.0 / sample_rate_hz;
    hp_alpha = (float) (rc / (rc + dt));	// For y[i] = alpha * (y[i-1] + x[i] - x[i-1])

    hp_x_prev = 0.0f;
    hp_y_prev = 0.0f;
}

float hipass_filter_process (float input_sample)
{
    if (hp_alpha == 1.0f)
	return input_sample;	// Filter disabled or pass-through

    float output_sample = hp_alpha * (hp_y_prev + input_sample - hp_x_prev);
    hp_y_prev = output_sample;
    hp_x_prev = input_sample;
    return output_sample;
}

void apply_hipass_filter_to_buffer (float *buffer, int num_samples, float cutoff_freq_hz, float sample_rate_hz)
{
    hipass_filter_setup (cutoff_freq_hz, sample_rate_hz);
    int i;
    for (i = 0; i < num_samples; i++)
    {
	buffer[i] = hipass_filter_process (buffer[i]);
    }
}

float returnMedianValue (float a, float b, float c, float d, float e)
{
    float sorted_arr[5];
    sorted_arr[0] = a;
    sorted_arr[1] = b;
    sorted_arr[2] = c;
    sorted_arr[3] = d;
    sorted_arr[4] = e;
    int i, j;
    float key;
    for (i = 1; i < 5; i++)
    {
	key = sorted_arr[i];
	j = i - 1;
	while (j >= 0 && sorted_arr[j] > key)
	{
	    sorted_arr[j + 1] = sorted_arr[j];
	    j = j - 1;
	}
	sorted_arr[j + 1] = key;
    }
    return sorted_arr[2];	// Median is the middle element (index 2 for 5 elements)
}

void save_all (void)
{
    save_encoder_state ();
    save_decoder_state ();
}

void restore_all (void)
{
    restore_encoder_state ();
    restore_decoder_state ();
}


void save_encoder_state (void)
{
    save_encoder_step_size = encoder_step_size;
    save_encoder_history = encoder_history;
}

void save_decoder_state (void)
{
    save_decoder_step_size = decoder_step_size;
    save_decoder_history = decoder_history;
}

void restore_encoder_state (void)
{
    encoder_step_size = save_encoder_step_size;
    encoder_history = save_encoder_history;
}

void restore_decoder_state (void)
{
    decoder_step_size = save_decoder_step_size;
    decoder_history = save_decoder_history;
}


void usage (char *programname)
{
    fprintf (stderr, "%s %s %s\n", PROGNAME, __DATE__, __TIME__);
    fprintf (stderr, "Usage: %s -i INPUTFILE -o OUTFILE [-r RATE | -R]\n", programname);
    fprintf (stderr, "\n");
    fprintf (stderr, "    Options for specifying input and output format details.\n");
    fprintf (stderr, "       -i specifies the input file. (WAV, MP3, OGG, etc.)\n");
    fprintf (stderr, "       -o specifies the output file name.\n");
    fprintf (stderr, "       -r specifies the output bitrate in samples/second. (Default: %ld)\n", TARGETRATE);
    fprintf (stderr, "       -R uses the input sample bitrate as the encoding rate.\n");
    fprintf (stderr, "\n");
}
