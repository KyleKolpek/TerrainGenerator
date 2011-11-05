///////////////////////////////////////////////////////////////////////////////
// ImprovedPerlinNoise
// Generates noise as described by Perlin at:
// http://mrl.nyu.edu/~perlin/paper445.pdf
// http://www.noisemachine.com/talk1/15.html
// http://cs.nyu.edu/~perlin/noise/
// This is not meant to be an efficient implementation. It is meant as a
// relatively verbose tool for pre-rendering noise. It will produce
// different results each time it is run.
///////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <bitset>
#include "PNGTools.h"

using namespace std;

void generateGradients(int dimensions);
int getGradientIndexAtCoords(
    int *const permutations,
    int permutationCount,
    int *const coords,
    int dimensions,
    int gradientCount);
double spline(
    double t);
double dot(
    double *a,
    double *b,
    double dimensions);
double interpolate(
    double state,
    double start,
    double end);
double noise(
    double *coords,
    int dimensions,
    int *const permutations,
    int permutationCount);

char **gVectors;
int edgeCount;

int main(
    int argc,
    char* argv[])
{
    if(argc < 6)
    {
        cerr << "ImprovedPerlinNoise.exe permutationCount bitDepth outFile "
                "gridSize x y [z w..]"<< endl;
        return 1;
    }
    int dimensions = argc-5;
    int permutationCount = strtol(argv[1],NULL,0); // May need to change this
    int bitDepth = strtol(argv[2],NULL,0);
    int gridSize = strtol(argv[4], NULL, 0);

    // Check for reasonable input
    if(bitDepth != 8 || bitDepth != 16)
    {
        cerr << "Error: bitDepth must be 8 or 16" << endl;
        return 1;
    }
    if(permutationCount <= 0)
    {
        cerr << "Error: permutationCount must be >= 1" << endl;
        return 1;
    }
    if(gridSize <= 0)
    {
        cerr << "Error: gridSize must be >= 1" << endl;
        return 1;
    }


    int *imageSize = new int[dimensions];
    for(int i=0; i<dimensions; i++)
    {
        imageSize[i] = strtol(argv[5+i], NULL, 0);
     
        if(imageSize[i] <= 0)
        {
            cerr << "Error: dimension " << i << " must be >= 1" << endl;
            return 1;
        }

    }

    // Calculate the number of edges in an n-dimensional hypercube
    edgeCount  = dimensions*pow(2.0, dimensions-1);

    // Initialize gradients
    generateGradients(dimensions);

    ///////////////////////////////////////////////////////////////////////////
    // Initialize permutations
    ///////////////////////////////////////////////////////////////////////////
    int *permutations = new int[permutationCount];
    srand(time(NULL));
    for(int i=0; i<permutationCount; i++)
    {
        permutations[i] = rand();
    }
    
    ///////////////////////////////////////////////////////////////////////////
    // Get Normalized Data
    ///////////////////////////////////////////////////////////////////////////

    // Assume 2D for now
    if(dimensions!=2)
    {
        return 1;
    }

    png_uint_16 **normalizedMap16;
    png_byte **normalizedMap8;
    
    // Create rowPtrs that wil be used for writing the image to a file later
    png_byte **rowPtrs = new png_byte*[imageSize[1]];

    if( bitDepth == 8 )
    {
        normalizedMap8 = new png_byte*[imageSize[1]];
    }
    else if( bitDepth == 16 )
    {
        normalizedMap16 = new png_uint_16*[imageSize[1]];
    }

    for (int r=0; r < imageSize[1]; r++)
    {
        if( bitDepth == 8 )
        {
            normalizedMap8[r] = new png_byte[imageSize[0]];
            rowPtrs[r] = (png_bytep)(normalizedMap8[r]);
        }
        else if( bitDepth == 16 )
        {
            normalizedMap16[r] = new png_uint_16[imageSize[0]];
            rowPtrs[r] = (png_bytep)(normalizedMap16[r]);
        }
    }
    
    // Calculate the maximum/minimum possible results of the dot product
    // operation and use it as the range when normalizing
    // This could be tuned since we are working with a restricted set of
    // vectors
    double limit = sqrt((double)dimensions*(dimensions-1));

    // Store normalized values
    for (int r=0; r < imageSize[1]; r++)
    {
        for (int c=0; c < imageSize[0]; c++)
        {
            // Translate pixel coords to grid coords
            double coords[2];
            coords[0] = (float)gridSize/imageSize[0] * c;
            coords[1] = (float)gridSize/imageSize[1] * r;

            if(bitDepth == 8)
            {
                normalizedMap8[r][c] = PNGTools::normalizeTo8Bit(
                                                     (noise(coords,
                                                           dimensions,
                                                           permutations,
                                                           permutationCount)
                                                      + limit) / (2*limit));
            }
            else if(bitDepth == 16)
            {
                normalizedMap16[r][c] = PNGTools::normalizeTo16Bit(
                                                      (noise(coords,
                                                           dimensions,
                                                           permutations,
                                                           permutationCount)
                                                      + limit) / (2*limit));
            }
        }
    }
    
    // Write to the file
    PNGTools::writePNGFile(argv[3],
                           rowPtrs,
                           imageSize[0],
                           imageSize[1],
                           bitDepth);

    ///////////////////////////////////////////////////////////////////////////
    // Cleanup
    ///////////////////////////////////////////////////////////////////////////

    for(int i=0; i < edgeCount; i++)
    {
        delete[] gVectors[i];
    }
    delete[] gVectors;

    delete[] permutations;

    for (int r=0; r<imageSize[1]; r++)
    {
        if(bitDepth == 8)
        {
            delete[] normalizedMap8[r];
        }
        else if(bitDepth == 16)
        {
            delete[] normalizedMap16[r];
        }
    }
    delete[] rowPtrs;

    delete[] imageSize;
    return 0;
}

///////////////////////////////////////////////////////////////////////////////
// generateGradients()
// Generate gradient vectors to the center of the edges of
// an n-dimensional cube
// Stores them in global gVectors for now
///////////////////////////////////////////////////////////////////////////////
void generateGradients(int dimensions)
{
    // An n-dimensional cube has n*2^(n-1) edges
    gVectors = new char*[edgeCount];
    for(int i=0; i < edgeCount; i++)
    {
        gVectors[i] = new char[dimensions];
    }

    // Generate permutations of (+-1,+-1,..,+-1) across n dimensions
    // One of the values should be a zero in order to point to the center of an
    // edge
    int index=0;
    for(int i=0; i < dimensions; i++)
    {
        for(int j=0; j < pow(2.0, (dimensions-1)); j++)
        {
            // Get our permutation information from j since we only need
            // binary values
            int permutation = j;
            for(int k=0; k<dimensions; k++)
            {
                if(i == k)
                {
                    gVectors[index][k]=0;
                }
                else
                {
                    // If the one bit of permutation is true assign 1, else -1
                    gVectors[index][k] = (permutation & 1) ? 1 : -1;
                    permutation        = permutation >> 1;
                }
            }
            index++;
        }

    }
}

///////////////////////////////////////////////////////////////////////////////
// getGradientIndexAtCoords()
// Produces a random index in [0-gradientCount)
// Guarantees the same output, given the same input
// randoms must consist of gradientCount random items
///////////////////////////////////////////////////////////////////////////////
int getGradientIndexAtCoords(
    int *const permutations,
    int permutationCount,
    int *const coords,
    int dimensions,
    int gradientCount)
{
    // Only work if we are given two or greater dimensions
    if(dimensions < 2)
    {
        return -1;
    }

    // Get a hash value
    int value = 0;
    for(int i=0; i<dimensions; i++)
    {
        value = permutations[(value + coords[i]) % permutationCount];
    }
    return value % gradientCount;
}

///////////////////////////////////////////////////////////////////////////////
// spline()
// Produces a weight given a state t using the equation 6t^5-15t^4+10t^3
///////////////////////////////////////////////////////////////////////////////
double spline(
    double t)
{
    // We'll just copy Perlin's equation
    return t*t*t*(t*(t*6-15)+10);
}

///////////////////////////////////////////////////////////////////////////////
// dot()
// Computes the dot product of two n-dimensional vectors: a and b
// Could easily be optimized since we are only ever multiplying 0,1,-1
///////////////////////////////////////////////////////////////////////////////
double dot(
    char *a,
    double *b,
    double dimensions)
{
    double sum = 0;
    for(int i=0; i<dimensions; i++)
    {
        sum += a[i]*b[i];
    }
    return sum;
}

///////////////////////////////////////////////////////////////////////////////
// interpolate()
// Linear interpolation
///////////////////////////////////////////////////////////////////////////////
double interpolate(
    double state,
    double start,
    double end)
{
    return start + state * (end - start);
}

///////////////////////////////////////////////////////////////////////////////
// noise()
// Takes coordinates and computes the noise at that point
// Completely deterministic
///////////////////////////////////////////////////////////////////////////////
double noise(
    double *coords,
    int dimensions,
    int *const permutations,
    int permutationCount)
{
    // Process coordinates
    int *hCoords = new int[dimensions];
    double *rCoords = new double[dimensions],
           *sWeights = new double[dimensions];
    for(int i=0; i<dimensions; i++)
    {
        // Get coordinates of containing hypercube
        hCoords[i] = (int)floor(coords[i]);

        // Get coordinates relative to containing hypercube [0-1)
        rCoords[i] = coords[i] - hCoords[i];

        // Get the weight of the spline at the relative coords
        sWeights[i] = spline(rCoords[i]);
    }



    ///////////////////////////////////////////////////////////////////////////
    // Calculate linear funciton
    ///////////////////////////////////////////////////////////////////////////

    // For every corner calculate the dot product of the associated
    // gradient and the relative position vector
    int numCorners = (int)pow(2.0,dimensions);
    double *dotResults = new double[numCorners];
    for(int i=0; i<numCorners; i++)
    {
        // Get the appropriate coordinates by using i to generate a
        // permutation and modify rCoords and hCoords accordingly
        int perm = i;
        for(int j=0; j<dimensions; j++)
        {
            rCoords[j] -= (perm >> j) & 1;
            hCoords[j] += (perm >> j) & 1;
        }
        int index = getGradientIndexAtCoords(permutations,
                                             permutationCount,
                                             hCoords,
                                             dimensions,
                                             edgeCount);
        dotResults[i] = dot(gVectors[index],
                            rCoords,
                            dimensions);
        
        // Undo the modifications of rCoords and hCoords
        for(int j=0; j<dimensions; j++)
        {
            rCoords[j] += (perm >> j) & 1;
            hCoords[j] -= (perm >> j) & 1;
        }
    }
    
    ///////////////////////////////////////////////////////////////////////////
    // Perform interpolation
    ///////////////////////////////////////////////////////////////////////////
    for(int i=1; i<=dimensions; i++)
    {
        for(int j=0; j<numCorners; j+=(int)pow(2.0,i) )
        {
            dotResults[j] = interpolate(sWeights[i-1],
                                        dotResults[j],
                                        dotResults[j+i]);
        }
    }

    double rVal = dotResults[0];

    ///////////////////////////////////////////////////////////////////////////
    // Cleanup
    ///////////////////////////////////////////////////////////////////////////
    delete[] hCoords;
    delete[] rCoords;
    delete[] sWeights;
    delete[] dotResults;

    return rVal;
}
