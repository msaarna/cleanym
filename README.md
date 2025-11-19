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
CLEANYM look-ahead is entirely implemented in the encoder, with a stock YMZ decoder, allowing for usage in applications with YMZ specific hardware.

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

## Reference Implementations

For the C reference implementation of the decoder logic, please refer to the `decodesample_ym4` function in `cleanym_decoder.c`.

Below are optimized assembly implementations for decoding a single byte (2 samples) on ARM processors.

### Shared Data Section
```asm
.data
.align 2
    @ The decoder state variables.
    @ Initialize these to: step=127, history=0 before first use.
decoder_state:
    .short  127     @ step_size (offset 0)
    .short  0       @ history   (offset 2)

    @ Step adaptation table (taken from C source)
    @ 230, 230, 230, 230, 307, 409, 512, 614
step_table:
    .word   230, 230, 230, 230, 307, 409, 512, 614
```

### ARMv7 Assembly (ARM State)
*Optimized for performance on Cortex-A / Cortex-R.*

```asm
.text
.arm
.align 4
.global decode_byte_ym4_arm

@ R0: Input Byte, R1: Output Buffer (int16_t*)
decode_byte_ym4_arm:
    PUSH    {R4-R8, LR}
    LDR     R4, =decoder_state
    LDRSH   R2, [R4, #0]        @ R2 = step
    LDRSH   R3, [R4, #2]        @ R3 = history
    LDR     R8, =step_table     @ R8 = table base
    MOV     R7, #2              @ Loop counter (2 nibbles per byte)

loop_arm:
    AND     R5, R0, #0x0F       @ Extract 4-bit sample (nibble)
    LSR     R0, R0, #4          @ Shift input for next iteration

    @ --- Calculate Diff ---
    @ diff = ((1 + (delta << 1)) * step) >> 3
    AND     R6, R5, #7          @ R6 = delta
    ADD     R6, R6, R6          @ R6 = delta * 2
    ADD     R6, R6, #1          @ R6 = 1 + (delta * 2)
    MUL     R6, R6, R2          @ R6 = operand * step
    LSR     R6, R6, #3          @ R6 = diff

    @ --- Update History ---
    TST     R5, #8              @ Check sign bit (bit 3)
    SUBNE   R3, R3, R6          @ If Sign, history -= diff
    ADDEQ   R3, R3, R6          @ Else, history += diff
    SSAT    R3, #16, R3         @ Clamp History to 16-bit signed

    STRH    R3, [R1], #2        @ Store Output Sample

    @ --- Update Step Size ---
    @ nstep = (step_table[delta] * step) >> 8
    AND     R6, R5, #7          @ R6 = delta (index)
    LDR     R6, [R8, R6, LSL #2]@ Load multiplier from table
    MUL     R2, R2, R6
    LSR     R2, R2, #8          @ step >> 8

    @ --- Clamp Step Size ---
    @ Range: 127 - 24576
    CMP     R2, #127
    MOVLT   R2, #127
    MOVW    R6, #24576
    CMP     R2, R6
    MOVGT   R2, R6

    SUBS    R7, R7, #1
    BNE     loop_arm

    @ Save state back to memory
    STRH    R2, [R4, #0]
    STRH    R3, [R4, #2]
    POP     {R4-R8, PC}
```

### Thumb-2 Assembly
*Optimized for code density on Cortex-M.*

```asm
.text
.thumb
.syntax unified
.align 2
.global decode_byte_ym4_thumb

@ R0: Input Byte, R1: Output Buffer (int16_t*)
decode_byte_ym4_thumb:
    PUSH    {R4-R8, LR}
    LDR     R4, =decoder_state
    LDRSH   R2, [R4, #0]        @ R2 = step
    LDRSH   R3, [R4, #2]        @ R3 = history
    LDR     R8, =step_table     @ R8 = table base
    MOV     R7, #2              @ Loop counter

loop_thumb:
    AND     R5, R0, #0x0F       @ Extract nibble
    LSR     R0, R0, #4          @ Shift input

    @ --- Calculate Diff ---
    AND     R6, R5, #7          @ delta
    ADD     R6, R6, R6          @ delta * 2
    ADD     R6, R6, #1          @ 1 + (delta * 2)
    MUL     R6, R6, R2          @ operand * step
    LSR     R6, R6, #3          @ diff

    @ --- Update History ---
    TST     R5, #8              @ Sign bit
    ITE     NE
    SUBNE   R3, R3, R6
    ADDEQ   R3, R3, R6
    SSAT    R3, #16, R3         @ Clamp

    STRH    R3, [R1], #2        @ Store

    @ --- Update Step Size ---
    AND     R6, R5, #7          @ delta
    LDR     R6, [R8, R6, LSL #2]@ Load multiplier
    MUL     R2, R2, R6
    LSR     R2, R2, #8          @ step >> 8

    @ --- Clamp Step Size ---
    CMP     R2, #127
    IT      LT
    MOVLT   R2, #127
    MOVW    R6, #24576
    CMP     R2, R6
    IT      GT
    MOVGT   R2, R6

    SUBS    R7, R7, #1
    BNE     loop_thumb

    STRH    R2, [R4, #0]
    STRH    R3, [R4, #2]
    POP     {R4-R8, PC}
```
