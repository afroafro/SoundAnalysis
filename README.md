# SoundAnalysis

This is a test project to make Qt GUI software (C++) that opens/saves sound files (.wav) and analyzes sound structures.

## Prerequisites

### Dependencies

 - Qt (Qt 5.7.0 (MSVC 2013, 32 bit))
 - FFTW (fftw-3.3.5-dll32) (QWidgetWAVChart.cpp includes it)

Folders should be set up in the parent's parent directory.

### Installation

Maybe I should prepare CMake. 

```
but I am not quite sure how to do that properly
```

I am using windows Qt Creator (MinGW 32bit) to build.

## Usage: 

Build them and you will see a main window pops up. Open directory, then listitems show names of sound files (.wav) in the directory. You can click file name to show its raw data and spectrogram.

## License

Free software under the terms of the GNU General Public License.
