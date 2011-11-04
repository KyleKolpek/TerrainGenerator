#include "PNGTools.h"
#include <cstdint>

namespace PNGTools
{
    png_byte normalizeTo8Bit(double original)
    {
        if(original >= 1)
            original = .9999;
        if(original < 0)
            original = 0;
        return (png_byte)(original*((png_byte)(-1)));
    }

    png_uint_16 normalizeTo16Bit(double original)
    {
        if(original >= 1)
            original = .9999;
        if(original < 0)
            original = 0;
        return (png_uint_16)(original*((png_uint_16)(-1)));

    }
    
    void writePNGFile(char* file_name, png_bytep* rowPtrs, int columns, int rows, png_byte bit_depth)
    {
        png_byte color_type=0;
        png_structp png_ptr;
        png_infop info_ptr;

        /* create file */
        FILE *fp = fopen(file_name, "wb");
        if (!fp)
            return;


        /* initialize stuff */
        png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
        //turn of compression
        png_set_compression_level(png_ptr,0);

        if (!png_ptr)
            return;

        info_ptr = png_create_info_struct(png_ptr);
        if (!info_ptr)
            return;

        if (setjmp(png_jmpbuf(png_ptr)))
            return;

        png_init_io(png_ptr, fp);

    
        //test
        png_set_compression_level(png_ptr, 0);
        png_set_filter(png_ptr, 0, PNG_FILTER_NONE);
        png_set_rows(png_ptr,info_ptr,rowPtrs);

        png_set_IHDR(png_ptr, info_ptr, columns, rows,
                 bit_depth, color_type, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
    
        //png_write_info(png_ptr, info_ptr);


        /* write bytes */
        //png_write_image(png_ptr, rowPtrs);


        // Check if we need to swap endian
        unsigned int transforms = 0;
        if(!isBigEndian())
        {
            transforms = PNG_TRANSFORM_SWAP_ENDIAN;
        }
        png_write_png(png_ptr,info_ptr,transforms,NULL);

        /* end write */
        //png_write_end(png_ptr, NULL);
    
        fclose(fp);
    }

    bool isBigEndian()
    {
        union {
            uint32_t i;
            char c[4];
        } x = {0x01020304};

        return x.c[0] == 1; 
    }
}