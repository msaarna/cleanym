# CLEANYM ADPCM Encoder/Decoder

CLEANYM is a project that implements a 4-bit YMZ ADPCM encoder and decoder. The encoder features an advanced look-ahead optimization that significantly improves audio quality compared to standard ADPCM algorithms and the stock YM ADPCM, especially for devices with low CPU capabilities.

The output is a raw stream of 4-bit ADPCM samples, without any specific header or container format, making it ideal for direct integration into embedded systems or applications requiring a minimal footprint.

## Features

*   **4-bit YMZ ADPCM Encoding:** Efficient compression suitable for low-bandwidth audio.
*   **Look-Ahead Optimization:** The encoder analyzes future samples to make more intelligent encoding decisions, significantly reducing perceived distortion and improving overall sound fidelity.
*   **Raw Sample Stream Output:** Generates a direct stream of ADPCM nibbles, packed into bytes, eliminating overhead for header parsing.
*   **High-Quality Decoding:** A robust decoder to convert the 4-bit ADPCM stream back into PCM audio.
*   **Resampling Support:** The encoder can resample input audio to a target rate, or preserve the original sample rate.
*   **Psychoacoustic Error Metric:** The look-ahead system uses a weighted psychoacoustic model to prioritize errors that are more noticeable to the human ear, focusing on absolute error and volatility (changes in error).

## Why CLEANYM?

Traditional ADPCM encoders often make decisions based only on the current sample, leading to "greedy" choices that can propagate errors or introduce noticeable artifacts. CLEANYM's look-ahead mechanism allows the encoder to anticipate future signal changes and choose the ADPCM code that minimizes distortion over a short window, resulting in a cleaner reproduction of the original audio. This is particularly beneficial for applications where every bit of quality matters within tight constraints.
CLEANYM look-ahead is entirely implemented in the encoder, with a stock YMZ decoder, allwoing for usage in applications with YMZ specific hardware.

## Build Instructions

To build CLEANYM, you will need a C compiler (like GCC or Clang) and the `libsndfile` and `libsamplerate` development libraries.

On Debian/Ubuntu:

```bash
sudo apt-get update
sudo apt-get install build-essential libsndfile1-dev libsamplerate0-dev
```

Then, compile the encoder and decoder:

```bash
gcc cleanym_encoder.c -o cleanym_encoder -lsndfile -lsamplerate -lm
gcc cleanym_decoder.c -o cleanym_decoder -lsndfile -lm
```

## Usage

### Encoder (`cleanym_encoder`)

Encodes a WAV file (or other `libsndfile`-supported formats) into a raw 4-bit YMZ ADPCM stream.

```bash
./cleanym_encoder -i <input_audio_file> -o <output_adpcm_file> [-r <rate> | -R]
```

**Options:**

*   `-i <input_audio_file>`: Specifies the input audio file (e.g., `input.wav`). **Required.**
*   `-o <output_adpcm_file>`: Specifies the output raw ADPCM file (e.g., `output.adpcm`). **Required.**
*   `-r <rate>`: Specifies the output sample rate in samples/second. Default is `32160`.
*   `-R`: Uses the input file's sample rate as the encoding rate, overriding `-r`.

**Example:**

```bash
./cleanym_encoder -i my_song.wav -o my_song.adpcm -r 16000
```
This will encode `my_song.wav` to a 16kHz 4-bit ADPCM stream and save it as `my_song.adpcm`.

### Decoder (`cleanym_decoder`)

Decodes a raw 4-bit YMZ ADPCM stream back into a WAV file.

```bash
./cleanym_decoder -i <input_adpcm_file> -o <output_wav_file> [-r <rate>]
```

**Options:**

*   `-i <input_adpcm_file>`: Specifies the input raw ADPCM file. **Required.**
*   `-o <output_wav_file>`: Specifies the output WAV file. **Required.**
*   `-r <rate>`: Specifies the sample rate of the output WAV file. This **must match** the rate used during encoding. Default is `14000`.

**Example:**

```bash
./cleanym_decoder -i my_song.adpcm -o decoded_my_song.wav -r 16000
```
This will decode `my_song.adpcm` (assuming it was encoded at 16kHz) and save it as `decoded_my_song.wav`.

## ADPCM Stream Format

The output ADPCM file is a simple byte stream where each byte contains two 4-bit ADPCM samples (nibbles). The lower nibble (bits 0-3) corresponds to the first sample, and the upper nibble (bits 4-7) corresponds to the second sample.

For example, if the encoder outputs `0xAB` for a byte, this represents two ADPCM samples: `0xB` and `0xA`. The decoder processes them in order: `0xB` then `0xA`.

If the total number of samples is odd, the last byte will contain one valid ADPCM nibble in the lower half, and the upper half will be zero.

