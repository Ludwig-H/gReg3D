#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

#include "Config.h"
#include "GDelaunay.h"
#include "GDelHost.h"
#include "Geometry.h"

// Scale a coordinate to the grid used by gReg3D
static RealType scalePoint( RealType gridWidth, float minVal, float maxVal, RealType inVal )
{
    inVal = inVal - minVal; // translate
    const float rangeVal = maxVal - minVal;
    inVal = ( gridWidth - 3.0f ) * inVal / rangeVal; // scale
    inVal += 1.0f; // shift to [1, gridWidth-2]
    return inVal;
}

int main( int argc, char** argv )
{
    if ( argc != 3 )
    {
        std::cerr << "Usage: EdgesWeightedDelaunay3D <input.xyzw> <output.txt>\n";
        return 1;
    }

    const std::string inFilename  = argv[1];
    const std::string outFilename = argv[2];

    std::ifstream inFile( inFilename.c_str() );
    if ( !inFile )
    {
        std::cerr << "Cannot open input file: " << inFilename << "\n";
        return 1;
    }

    struct RawPoint { float x, y, z, w; };
    std::vector< RawPoint > rawPts;
    rawPts.reserve( 1024 );

    float minVal = std::numeric_limits<float>::max();
    float maxVal = std::numeric_limits<float>::lowest();
    float x, y, z, w;

    while ( inFile >> x >> y >> z >> w )
    {
        rawPts.push_back( { x, y, z, w } );
        if ( x < minVal ) minVal = x;
        if ( y < minVal ) minVal = y;
        if ( z < minVal ) minVal = z;
        if ( x > maxVal ) maxVal = x;
        if ( y > maxVal ) maxVal = y;
        if ( z > maxVal ) maxVal = z;
    }
    inFile.close();

    Point3HVec pointVec;
    WeightHVec weightVec;
    pointVec.reserve( rawPts.size() );
    weightVec.reserve( rawPts.size() );

    const RealType gridSize = 512.0f;
    for ( const RawPoint& rp : rawPts )
    {
        Point3 p;
        p._p[0] = scalePoint( gridSize, minVal, maxVal, rp.x );
        p._p[1] = scalePoint( gridSize, minVal, maxVal, rp.y );
        p._p[2] = scalePoint( gridSize, minVal, maxVal, rp.z );
        pointVec.push_back( p );
        weightVec.push_back( rp.w );
    }

    Config config;
    config._run        = 0;
    config._runNum     = 1;
    config._gridSize   = (int)gridSize;
    config._pointNum   = static_cast<int>( pointVec.size() );
    config._dist       = UniformDistribution;
    config._facetMax   = 12000000; // default upper bound
    config._weightMax  = 1;
    config._logVerbose = false;
    config._logStats   = false;
    config._logTiming  = false;
    config._doCheck    = false;
    config._inFile     = true;
    config._inFilename = inFilename;

    gdelInit( config, pointVec, weightVec );

    double t1, t2, t3, t4;
    gdelCompute( t1, t2, t3, t4 );

    const TetraHVec& tetraVec = getHostTetra();

    SegmentHVec segVec;
    segVec.reserve( tetraVec.size() * 6 );
    Segment segArr[6];
    for ( const Tetrahedron& tet : tetraVec )
    {
        tet.getSegments( segArr );
        std::copy( segArr, segArr + 6, std::back_inserter( segVec ) );
    }

    std::sort( segVec.begin(), segVec.end() );
    segVec.erase( std::unique( segVec.begin(), segVec.end() ), segVec.end() );

    std::ofstream outFile( outFilename.c_str() );
    if ( !outFile )
    {
        std::cerr << "Cannot open output file: " << outFilename << "\n";
        gdelDeInit();
        return 1;
    }

    for ( const Segment& s : segVec )
    {
        outFile << s._v[0] << " " << s._v[1] << "\n";
    }
    outFile.close();

    gdelDeInit();
    return 0;
}

