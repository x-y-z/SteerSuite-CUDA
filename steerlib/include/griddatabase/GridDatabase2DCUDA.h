#ifndef __GRIDDATABASE2DCUDA_H__
#define __GRIDDATABASE2DCUDA_H__

#include "cuda.h"
#include "cuda_runtime.h"

#include "testcaseio/ObstacleInitialConditions.h"
#include <cstdio>
#include "util/Geometry.h"
#include "cutil_math.h"


#define STREAM_SIZE 512

//=============cuda===========
// Enable this for error checking
#define CUDA_CHECK_ERROR

#define CudaSafeCall( err )     __cudaSafeCall( err, __FILE__, __LINE__ )
#define CudaCheckError()        __cudaCheckError( __FILE__, __LINE__ )

inline void __cudaSafeCall( cudaError err, const char *file, const int line )
{
#ifdef CUDA_CHECK_ERROR

#pragma warning( push )
#pragma warning( disable: 4127 ) // Prevent warning on do-while(0);

    do
    {
        if ( cudaSuccess != err )
        {
            fprintf( stderr, "cudaSafeCall() failed at %s:%i : %s\n",
                     file, line, cudaGetErrorString( err ) );
            exit( -1 );
        }
    } while ( 0 );

#pragma warning( pop )

#endif  // CUDA_CHECK_ERROR

    return;
}

inline void __cudaCheckError( const char *file, const int line )
{
#ifdef CUDA_CHECK_ERROR

#pragma warning( push )
#pragma warning( disable: 4127 ) // Prevent warning on do-while(0);

    do
    {
        cudaError_t err = cudaGetLastError();
        if ( cudaSuccess != err )
        {
            fprintf( stderr, "cudaCheckError() failed at %s:%i : %s.\n",
                     file, line, cudaGetErrorString( err ) );
            exit( -1 );
        }

        // More careful checking. However, this will affect performance.
        // Comment if not needed.
        err = cudaThreadSynchronize();
        if( cudaSuccess != err )
        {
            fprintf( stderr, "cudaCheckError() with sync failed at %s:%i : %s.\n",
                     file, line, cudaGetErrorString( err ) );
            exit( -1 );
        }
    } while ( 0 );

#pragma warning( pop )

#endif // CUDA_CHECK_ERROR

    return;
}

//=============cuda end=========

typedef struct AABox{
	float xmin, xmax;
	float ymin, ymax;
	float zmin, zmax;

	AABox &operator= (const Util::AxisAlignedBox &aBox)
	{
		//AABox mBox;
		xmax = aBox.xmax;
		xmin = aBox.xmin;
		ymax = aBox.ymax;
		ymin = aBox.ymin;
		zmax = aBox.zmax;
		zmin = aBox.zmin;

		return *this;
	}

	AABox &operator= (const SteerLib::ObstacleInitialConditions &aBox)
	{
		//AABox mBox;
		xmax = aBox.xmax;
		xmin = aBox.xmin;
		ymax = aBox.ymax;
		ymin = aBox.ymin;
		zmax = aBox.zmax;
		zmin = aBox.zmin;

		return *this;
	}
}AABox;


typedef struct cuda_agent{
	bool _enabled;
	float3 _position;
	float3 _velocity;
	float3 _forward;
	float _radius;
	float3 _goalQueue[20];
	int _goalQSize;
	int _curGoal;
	int _usedGoal;
	AABox _oldBounds;
	AABox _newBounds;
}cuda_agent;




typedef struct cuda_obstacle{
	AABox _bounds;
	bool _blocksLineOfSight;
} cuda_obstacle;

typedef struct cuda_item{
	int type;//0: agent, 1: obstacle
	cuda_agent _agent;
	cuda_obstacle _obstacle;
} cuda_item;



void launch_updateAICUDA(cuda_item *cudaItems, cuda_item *hostItems, float currentSimulationTime, float simulatonDt, unsigned int currentFrameNumber, 
	        int agentNum, int obstacleNum, int &numDisabledAgents, int streamNum, cudaStream_t *stList);

#endif