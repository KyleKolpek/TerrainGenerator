#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <map>
#include <cfloat>
#include "png.h"
#include "PNGTools.h"

using namespace std;

typedef struct Point
{
    int r;
    int c;
    
    // Constructor
    Point(int r, int c)
    {
        this->r = r;
        this->c = c;
    }
}Point;

double getRelativeDistance( int x1, int y1, int x2, int y2 );
double getRealDistance( int x1, int y1, int x2, int y2 );

int main( int argc, char *argv[] )
{
    ///////////////////////////////////////////////////////////////////////////
    // Get input
    ///////////////////////////////////////////////////////////////////////////

    if(argc < 6)
    {
        cerr << "VoronoiDiagram.exe width height bit_depth point_count outfile c1 [c2 c3..]" << endl;
        return 1;
    }
    
    int columns = strtol(argv[1],NULL,0);
    int rows = strtol(argv[2],NULL,0);
    int bitDepth = strtol(argv[3],NULL,0);
    int pointCount = strtol(argv[4],NULL,0);
    int cCount = argc - 6;
    double *cVals = new double[cCount];

    for( int i=0; i<cCount; i++ )
    {
        cVals[i] = strtod(argv[6+i],NULL);
    }
    
    // Roughly check inputs
    if( columns < 1 )
    {
        cerr << "Error: width must be >= 1" << endl;
        return 1;
    }
    if( columns < 1 )
    {
        cerr << "Error: height must be >= 1" << endl;
        return 1;
    }
    if( pointCount < 1 )
    {
        cerr << "Error: point_count must be >= 1" << endl;
        return 1;
    } 
    if( cCount > pointCount)
    {
        cerr << "Error: number of c paramaters cannot exceed point_count" << endl;
    }
    ///////////////////////////////////////////////////////////////////////////
    // Generate Diagram
    ///////////////////////////////////////////////////////////////////////////

    // Choose random site points bounded by [(0,0),(width,height))
    srand(time(NULL));
    Point **points = new Point*[pointCount];
    for( int i=0; i<pointCount; i++ )
    {
        points[i] = new Point( rand()%rows, rand()%columns );
    }

    // Track largest and smallest values [Sum(cX*dX)] for normalization
    double largestVal = -FLT_MAX;
    double smallestVal = FLT_MAX;

    // Generate values at every free point
    double **image = new double*[rows];
    for( int i=0; i<rows; i++ )
    {
        image[i] = new double[columns];
    }

    // For every free point (pixel)
    for( int r=0; r<rows; r++ )
    {
        for( int c=0; c<columns; c++ )
        {
            // Find the closest points
            multimap<double,Point*> checkedPoints;
            multimap<double,Point*>::iterator it;
            for( int i=0; i<pointCount; i++ )
            {
                double relativeDist = getRelativeDistance(
                    points[i]->r, points[i]->c, r, c );
                checkedPoints.insert( pair<double,Point*>(
                    relativeDist, points[i] ) );
            }
            it = checkedPoints.begin();

            // Associate a value with the pixel
            image[r][c] = 0.0;
            for( int i=0; i<cCount; i++ )
            {
                double realDist = getRealDistance(
                    (*it).second->r, (*it).second->c, r, c );
                image[r][c] += cVals[i] * realDist;
                advance(it,1);
            }

            // Track min and max
            largestVal = image[r][c] > largestVal ?
                image[r][c] : largestVal;
            smallestVal = image[r][c] < smallestVal ?
                image[r][c] : smallestVal;
        }
    }
    
    ///////////////////////////////////////////////////////////////////////////
    // Normalize Data
    ///////////////////////////////////////////////////////////////////////////

    png_uint_16 **normalizedMap16;
    png_byte **normalizedMap8;

    // Normalize to 0-.999
    for( int r=0; r<rows; r++ )
    {
        for( int c=0; c<columns; c++ )
        {
            image[r][c] = image[r][c] - smallestVal;
            image[r][c] = image[r][c] / (largestVal - smallestVal);
        }
    }

    // Create rowPtrs that wil be used for writing the image to a file later
    png_byte **rowPtrs = new png_byte*[rows];

    // Normalize to 0-2^bitDepth-1
    if( bitDepth == 8 )
    {
        normalizedMap8 = new png_byte*[rows];
    }
    else if( bitDepth == 16 )
    {
        normalizedMap16 = new png_uint_16*[rows];
    }

    // Allocate memory for normalized data and point rowPtrs to it
    for (int r=0; r<rows; r++)
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
                normalizedMap8[r][c] = PNGTools::normalizeTo8Bit(image[r][c]);
            }
            else if(bitDepth == 16)
            {
                normalizedMap16[r][c] = PNGTools::normalizeTo16Bit(image[r][c]);
            }
        }
    }
    
    ///////////////////////////////////////////////////////////////////////////
    // Write out data
    ///////////////////////////////////////////////////////////////////////////
    
    PNGTools::writePNGFile(argv[5],rowPtrs,columns,rows,bitDepth);
    
    ///////////////////////////////////////////////////////////////////////////
    // Cleanup
    ///////////////////////////////////////////////////////////////////////////

    for( int i=0; i<pointCount; i++ )
    {
        delete points[i];
    }
    for (int r=0; r<rows; r++)
    {
        delete[] image[r];
        if(bitDepth == 8)
        {
            delete[] normalizedMap8[r];
        }
        else if(bitDepth == 16)
        {
            delete[] normalizedMap16[r];
        }
    }
    if(bitDepth == 8)
    {
        delete[] normalizedMap8;
    }
    else if(bitDepth == 16)
    {
        delete[] normalizedMap16;
    }
    delete[] image;
    delete[] points;
    delete[] cVals;
    delete[] rowPtrs;
    
    return 0;
}

// Useful for comparing distances between different points
// Doesn't use sqrt
double getRelativeDistance( int x1, int y1, int x2, int y2 )
{
    return (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);
}

// Used for getting the actual distance between two points
// More expensive than getRelativeDistance because of sqrt()
double getRealDistance( int x1, int y1, int x2, int y2 )
{
    // Don't use sqrt(getRelativeDistance(..)) because it may change
    return sqrt((double)((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)));
}