#include "griddatabase/GridDatabase2DCUDA.h"

#define BLOCKSIZE 16

#define MAX_FORCE_MAGNITUDE 3.0f
#define MAX_SPEED 1.3f
#define AGENT_MASS 1.0f

//cuda update agent kernel
__global__ void updateAI_kernel(cuda_item *cudaItems, float currentSimulationTime, float dt, unsigned int currentFrameNumbers,
								int agentNum, int obstacleNum, int *disabledAgents, int streamIdx)
{
	int x = blockIdx.x * blockDim.x + threadIdx.x;// + streamIdx*STREAM_SIZE;
	
	if (x >= agentNum)
		return;

	if (!cudaItems[x]._agent._enabled)
	{
		cudaItems[x].type = -1;
		disabledAgents[x] = 1;
		return;
	}
	//printf("thread: %d\t", x);

	int curGoal = cudaItems[x]._agent._curGoal;
	float radius = cudaItems[x]._agent._radius;
	float3 position = cudaItems[x]._agent._position;
	float3 vectorToGoal = cudaItems[x]._agent._goalQueue[curGoal] - cudaItems[x]._agent._position;

	
	cudaItems[x]._agent._usedGoal = 0;

	// it is up to the agent to decide what it means to have "accomplished" or "completed" a goal.
	// for the simple AI, if the agent's distance to its goal is less than its radius, then the agent has reached the goal.
	if (dot(vectorToGoal, vectorToGoal) < radius * radius) {
		cudaItems[x]._agent._curGoal++;
		cudaItems[x]._agent._usedGoal++;
		if (cudaItems[x]._agent._curGoal != cudaItems[x]._agent._goalQSize) {
			// in this case, there are still more goals, so start steering to the next goal.
			vectorToGoal = cudaItems[x]._agent._goalQueue[cudaItems[x]._agent._curGoal] - cudaItems[x]._agent._position;
		}
		else {
			// in this case, there are no more goals, so disable the agent and remove it from the spatial database.
			AABox bounds = {position.x-radius, position.x+radius, 0.0f, 0.0f, position.z-radius, position.z+radius};
			cudaItems[x]._agent._newBounds = bounds;
			cudaItems[x]._agent._enabled = false;
			disabledAgents[x] = 1;
			return;
		}
	}
	
	float3 clippedForce = clamp(vectorToGoal, MAX_FORCE_MAGNITUDE);
	float3 acceleration = (clippedForce / AGENT_MASS);
	cudaItems[x]._agent._velocity += (dt*acceleration);
	cudaItems[x]._agent._velocity = clamp(cudaItems[x]._agent._velocity, MAX_SPEED);  // clamp _velocity to the max speed
	float3 newPosition = cudaItems[x]._agent._position + (dt*cudaItems[x]._agent._velocity);

	// For this simple agent, we just make the orientation point along the agent's current velocity.
	if (dot(cudaItems[x]._agent._velocity,cudaItems[x]._agent._velocity) != 0.0f) {
		cudaItems[x]._agent._forward = normalize(cudaItems[x]._agent._velocity);
	}

	// update the database with the new agent's setup
	AABox oldBounds = {cudaItems[x]._agent._position.x - cudaItems[x]._agent._radius, 
					   cudaItems[x]._agent._position.x + cudaItems[x]._agent._radius, 
					   0.0f, 0.0f, 
					   cudaItems[x]._agent._position.z - cudaItems[x]._agent._radius, 
					   cudaItems[x]._agent._position.z + cudaItems[x]._agent._radius};
	AABox newBounds = {newPosition.x - cudaItems[x]._agent._radius, 
		               newPosition.x + cudaItems[x]._agent._radius, 
					   0.0f, 0.0f, 
					   newPosition.z - cudaItems[x]._agent._radius, 
					   newPosition.z + cudaItems[x]._agent._radius};

	cudaItems[x]._agent._oldBounds = oldBounds;
	cudaItems[x]._agent._newBounds = newBounds;

	cudaItems[x]._agent._position = newPosition;

}

void launch_updateAICUDA(cuda_item *cudaItems, cuda_item *hostItems, float currentSimulationTime, float simulatonDt, unsigned int currentFrameNumber, 
	                     int agentNum, int obstacleNum, int &numDisabledAgents, int streamNum, cudaStream_t *stList)
{
	dim3 block(BLOCKSIZE*BLOCKSIZE);
	dim3 grid((agentNum)/(BLOCKSIZE*BLOCKSIZE) + 1);
	//tried stream, but not work well
	//dim3 grid((STREAM_SIZE)/(BLOCKSIZE*BLOCKSIZE));

	int *disAgents, *hostDisAgents;
	CudaSafeCall(cudaMalloc(&disAgents, sizeof(int)*agentNum));
	CudaSafeCall(cudaMemset(disAgents,0, sizeof(int)*agentNum));

	hostDisAgents = new int[agentNum];
	//int streamAgentNum = 1024;

	//for (int i = 0; i < streamNum; ++i)
	//{
		//int streamSize;
	updateAI_kernel<<<grid, block>>>(cudaItems, currentSimulationTime, simulatonDt, currentFrameNumber,
		                            agentNum, obstacleNum, disAgents, 0);

	//for stream, not used
	/*if (STREAM_SIZE*(i+1) <= agentNum)
		streamSize = STREAM_SIZE;
	else
		streamSize = agentNum % STREAM_SIZE;

	CudaSafeCall(cudaMemcpyAsync(hostItems + STREAM_SIZE*i,
								cudaItems + STREAM_SIZE*i,
								sizeof(cuda_item)*streamSize,
								cudaMemcpyDeviceToHost,
								stList[i]));*/
	//}

	cudaError_t res;// = cudaDeviceSynchronize();

	res = (cudaMemcpy(hostDisAgents, disAgents, sizeof(int)*agentNum, cudaMemcpyDeviceToHost));

	for (int i = 0; i < agentNum; ++i)
	{
		numDisabledAgents += hostDisAgents[i];
	}

	CudaSafeCall(cudaFree(disAgents));

}