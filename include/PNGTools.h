#ifndef PNGTOOLS_H
#define PNGTOOLS_H

#include "png.h"

namespace PNGTools
{
	png_byte normalizeTo8Bit(double original);
	png_uint_16 normalizeTo16Bit(double original);
	void writePNGFile(char* file_name, png_bytep* rowPtrs, int columns, int rows, png_byte bit_depth);
	bool isBigEndian();
}

#endif