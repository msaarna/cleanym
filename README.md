# CLEANYM ADPCM Encoder/Decoder

CLEANYM is a project that implements a 4-bit YMZ ADPCM encoder and decoder. The encoder features an advanced look-ahead optimization that significantly improves audio quality compared to standard ADPCM algorithms and the stock YM ADPCM, especially for devices with low CPU capabilities.

The output is a raw stream of 4-bit ADPCM samples, without any specific header or container format, making it ideal for direct integration into embedded systems or applications requiring a minimal footprint.

## Features

*   **4-bit YMZ ADPCM Encoding:** Efficient compression suitable for low-bandwidth audio.
*   **Look-Ahead Optimization:** The encoder analyzes future samples to make more intelligent encoding decisions, drastically reducing perceived distortion and improving overall sound fidelity.
*   **Raw Sample Stream Output:** Generates a direct stream of ADPCM nibbles, packed into bytes, eliminating overhead for header parsing.
*   **High-Quality Decoding:** A robust decoder to convert the 4-bit ADPCM stream back into PCM audio.
*   **Resampling Support:** The encoder can resample input audio to a target rate, or preserve the original sample rate.
*   **Psychoacoustic Error Metric:** The look-ahead system uses a weighted psychoacoustic model to prioritize errors that are more noticeable to the human ear, focusing on absolute error and volatility (changes in error).

## Why CLEANYM?

Traditional ADPCM encoders often make decisions based only on the current sample, leading to "greedy" choices that can propagate errors or introduce noticeable artifacts. CLEANYM's look-ahead mechanism allows the encoder to anticipate future signal changes and choose the ADPCM code that minimizes distortion over a short window, resulting in a much cleaner and more faithful reproduction of the original audio. This is particularly beneficial for applications where every bit of quality matters within tight constraints.

## Build Instructions

To build CLEANYM, you will need a C compiler (like GCC or Clang) and the `libsndfile` and `libsamplerate` development libraries.

On Debian/Ubuntu:

```bash
sudo apt-get update
sudo apt-get install build-essential libsndfile1-dev libsamplerate0-dev
