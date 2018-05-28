# straight_process
A tool to process files with STRAIGHT.

[STRAIGHT](http://www.wakayama-u.ac.jp/~kawahara/STRAIGHTadv/index_e.html) is a high fidelity vocoder created by Hideki Kawahara at Wakayama University, Japan. It's been used in numerous research projects to manipulate vocal properties of speech recordings.

Since the initial version (STRAIGHTV40_006b), Kawahara and his colleagues have developed many new versions. This `straight_process` tool offers a common interface to a number of these tools.

In addition, the `straight_process` function caches the generated files (both from analysis and synthesis) to speed-up processing in the context of an experiment.

__Note__: STRAIGHT is not open source and cannot be freely distributed here. `straight_process` assumes that you've obtained STRAIGHT, Tandem-Straight, StraightLib or WORLD by you own means. STRAIGHT and Tandem-Straight can be obtained from this [webpage](http://www.wakayama-u.ac.jp/~kawahara/STRAIGHTadv/index_e.html). StraightLib does not seem to be available anymore at the moment. WORLD is available on [github](https://github.com/mmorise/World). Once you have (at least some of) the required libraries, you have to

## Usage

```matlab
STRAIGHT_PROCESS - Cached STRAIGHT process of files
  [Y, FS] = STRAIGHT_PROCESS(FILENAME, DF0)
      Processes FILENAME with STRAIGHTV40_006b to apply DF0, which
      represents:
      - If DF0 is a number, it is treated as the distance in semitones from the
        geometric mean F0.
      - If DF0 is a string of the form '100.1Hz', it is treated as the new
        absolue mean F0 in Hertz.
      - If DF0 is a string starting with '+' or '-' and ending with 'Hz',
        it is treated as a linear offset in Hertz.

  [Y, FS] = STRAIGHT_PROCESS(FILENAME, DF0, DVTL)
      DVTL is specified as a shift in semitones. ATTENTION: this is VTL
      change, i.e. if DVTL is positive, formant frequencies are shifted
      down.

  [Y, FS] = STRAIGHT_PROCESS(FILENAME, DF0, DVTL, DURATION)
      DURATION specifies the duration modificator:
      - If DURATION is a number, it specifies the new duration in
        seconds.
      - If DURATION is a string of the form '*1.5', the duration is
        multiplied by the number following '*'.
      - If DURATION is a string of the form '+1.5s' or '-1.5s', the duration is
        extended by the given number of seconds.
      - If DURATION is empty ([] or '') or NaN, the duration is unchanged.

  [Y, FS] = STRAIGHT_PROCESS(..., METHOD)
      METHOD is a string that is one of the following:
      - 'straight': then STRAIGHTV40_006b is used (Matlab).
      - 'tandem': then baseTandemSTRAIGHTV009x is used (Matlab).
      - 'straightlib': then straightLib 1.1.4 is used (binary).
      - 'world': Morise's WORLD vocoder (Matlab).

  [Y, FS] = STRAIGHT_PROCESS(..., PARAMS)
      Where PARAMS is a struct() giving extra parameters:
      - cache_folder: the path where the cached files are stored
      - cache_format: 'flac' (default) or 'wav'
      - straight_path, tandem_path, straightlib_path: the paths to the
        various STRAIGHT libraries.
      - lower_f0_bound, upper_f0_bound: these are the min and max F0
        values used during the F0 extraction process.

  [Y, FS] = STRAIGHT_PROCESS(X, FS, ...)
      Same as above but passing the signal directly instead of a file
      name (no cache is created).

  Important note on clipping: The returned signals are scaled such that
  their RMS is identical to the original sound. That may mean that their
  peak amplitude is larger than 1, which would result in clipping. The
  cache files are stored in a way that prevents clipping. It is the
  user's job to make sure that clipping does not occur at playback.
```

## Example

The call below will shift the F0 one octave down, extend the VTL by 3.8
semitones and double the duration:

```matlab
[y, fs] = straight_process('filename.wav', -12, 3.8, '*2', 'straight');
```

If the same wavfile is processed again, the analysis step will be skipped, and
instead the analysis will be read from the cached mat file.

