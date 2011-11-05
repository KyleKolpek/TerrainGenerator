///////////////////////////////////////////////////////////////////////////////
// DiamondSquare
// Produces fractal noise following the diamond-square algorithm.
// Currently restricted to dimensions (n^2)+1 width/height.
///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <cstdlib>
#include <ctime>
#include "png.h"
#include "PNGTools.h"

using namespace std;

double SEED = .5;

void square(
    double** heightMap,
    int rows,
    int columns,
    double variance,
    int sideLength);
void diamond(
    double** heightMap,
    int rows,
    int columns,
    double variance,
    int sideLength);
bool isPointValid(
    int rows,
    int columns,
    int row,
    int column);

int main(int argc, char* argv[]) {

    if(argc != 6)
    {
        cerr << "DiamondSquare.exe width/height bit_depth initial_variance"
                "rate_of_change(!=0) outfile" << endl;
        //      "[top_left=.5 top_right=.5 bottom_left=.5 bottom_right=.5]";
        return 1;
    }
    // Get command line parameters
    int columns  = strtol(argv[1], NULL, 0);
    int rows     = strtol(argv[1], NULL, 0);
    int bitDepth = strtol(argv[2], NULL, 0);

    //the range (-VARIANCE -> +VARIANCE) for the average offset
    double variance     = strtod(argv[3], NULL);
    double rateOfChange = strtod(argv[4], NULL);

    // Check if width/height is a (power of 2) + 1
    if( ((rows-1) & (rows-2)) != 0)
    {
        cerr << "ERROR: width/height must be a power of 2 + 1" << endl;
        return 1;
    }


    // heightMap contains values between 0 and .999
    double** heightMap;

    // normalizedMaps are used to store values once they have been normalized
    // to the range [0 - 2^bitDepth)
    png_uint_16 **normalizedMap16;
    png_byte **normalizedMap8;
    
    // rowPtrs contains pointers to the rows in normalizedMap
    // It is needed for writing the PNGs
    png_byte **rowPtrs;

    heightMap = new double*[rows];
    for (int i = 0; i < rows; i++)
        heightMap[i] = new double[columns];

    // Seed corners
    srand(time(NULL));
    heightMap[0][0] = (float)rand() /RAND_MAX;
    heightMap[0][columns - 1] = (float)rand() /RAND_MAX;
    heightMap[rows - 1][0] = (float)rand() /RAND_MAX;
    heightMap[rows - 1][columns - 1] = (float)rand() /RAND_MAX;

    // Perform diamond square algorithm
    for(int sideLength=columns-1;
            sideLength>1;
            sideLength/=2, variance*=rateOfChange)
    {
        square(heightMap, rows, columns, variance, sideLength);
        diamond(heightMap, rows, columns, variance, sideLength);
    }
    
    
    // Create normalized map
    if( bitDepth == 8 )
    {
        normalizedMap8 = new png_byte*[rows];
    }
    else if( bitDepth == 16 )
    {
        normalizedMap16 = new png_uint_16*[rows];
    }

    // Create row ptr array
    rowPtrs = new png_byte*[rows];

    for (int r=0; r < rows; r++)
    {
        if( bitDepth == 8 )
        {
            normalizedMap8[r] = new png_byte[columns];
            rowPtrs[r] = (png_bytep)(normalizedMap8[r]);
        }
        else if( bitDepth == 16 )
        {
            normalizedMap16[r] = new png_uint_16[columns];
            rowPtrs[r] = (png_bytep)(normalizedMap16[r]);
        }
    }

    // Store normalized values
    for (int r=0; r < rows; r++)
    {
        for (int c=0; c < columns; c++)
        {
            if(bitDepth == 8)
            {
                normalizedMap8[r][c] = PNGTools::normalizeTo8Bit(
                                           heightMap[r][c]);
            }
            else if(bitDepth == 16)
            {
                normalizedMap16[r][c] = PNGTools::normalizeTo16Bit(
                                            heightMap[r][c]);
            }
        }
    }

    // Write to the file
    PNGTools::writePNGFile(argv[5], rowPtrs, columns, rows, bitDepth);

    // Cleanup
    for (int r=0; r<rows; r++)
    {
        delete[] heightMap[r];
        if(bitDepth == 8)
        {
            delete[] normalizedMap8[r];
        }
        else if(bitDepth == 16)
        {
            delete[] normalizedMap16[r];
        }
    }
    delete[] heightMap;
    if(bitDepth == 8)
    {
        delete[] normalizedMap8;
    }
    else if(bitDepth == 16)
    {
        delete[] normalizedMap16;
    }
    delete[] rowPtrs;
}

///////////////////////////////////////////////////////////////////////////////
// square()
// Performs the "square" operations which takes the values at four corners of
// a square as described by the two corners [r,c] and
// [r+sideLength,c+sidelength] in heightMap and stores the average of the
// values + varience*(rand(-1,1)) in heightMap at the center of the four
// coordinates.
///////////////////////////////////////////////////////////////////////////////
void square(
    double** heightMap,
    int rows,
    int columns,
    double variance,
    int sideLength)
{
    int halfSide = sideLength/2;

    for(int r=0; r < rows-1; r+=sideLength)
    {
        for(int c=0; c < columns-1; c+=sideLength)
        {
            double avg = heightMap[r][c]+
                heightMap[r+sideLength][c+sideLength]+
                heightMap[r+sideLength][c]+
                heightMap[r][c+sideLength];
            avg /= 4;

            //fix overflow
            heightMap[r+halfSide][c+halfSide] =
                avg + (double)rand() / RAND_MAX * 2.0 * variance - variance;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
// diamond()
// Performs the "square" operations which takes the values at four corners of
// a diamond as described by the four corners [r+sideLength/2,c],
// [r,c+sidelength/2],[r+sideLength,c+sidelength/2],
// and [r+sidelength/2,c+sidelength] in heightMap and stores the average of the
// values + varience*(rand(-1,1)) in heightMap at the center of the four
// coordinates.
///////////////////////////////////////////////////////////////////////////////
void diamond(
    double** heightMap,
    int rows,
    int columns,
    double variance,
    int sideLength)
{
    int halfSide = sideLength/2;
    for(int c=0; c < columns; c+=halfSide)
    {
        for(int r=(c+halfSide)%sideLength; r < rows; r+=sideLength)
        {
            int count = 0;
            double avg = 0;
            if(isPointValid(rows, columns, r, c-halfSide))
            {
                count++;
                avg += heightMap[r][c-halfSide];
            }
            if(isPointValid(rows, columns, r, c+halfSide))
            {
                count++;
                avg += heightMap[r][c+halfSide];
            }
            if(isPointValid(rows, columns, r+halfSide, c))
            {
                count++;
                avg += heightMap[r+halfSide][c];
            }
            if(isPointValid(rows, columns, r-halfSide, c))
            {
                count++;
                avg += heightMap[r-halfSide][c];
            }
            avg /= count;
            //fix overflow
            heightMap[r][c]=
                avg + (double)rand() / RAND_MAX * 2.0 * variance - variance;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
// isPointValid()
// Simple check to make sure (row,column) is between (0,0) and (rows,columns).
// Returns true if this is the case, and false otherwise.
///////////////////////////////////////////////////////////////////////////////
bool isPointValid(
    int rows,
    int columns,
    int row,
    int column)
{
    if(row<0 || column<0)
    {
        return false;
    }
    if(row>=rows || column>=columns)
    {
        return false;
    }
    return true;
}

