//
// Copyright (c) 2009-2010 Shawn Singh, Mubbasir Kapadia, Petros Faloutsos, Glenn Reinman
// See license.txt for complete license.
//

#ifndef __SIMPLE_AGENT__
#define __SIMPLE_AGENT__

/// @file SimpleAgent.h
/// @brief Declares the SimpleAgent class.

#include <queue>
#include "SteerLib.h"
#include "SimpleAgent.h"
#include "SimpleAIModule.h"

/**
 * @brief An example agent with very basic AI, that is part of the simpleAI plugin.
 *
 * This agent performs extremely simple AI using forces and Euler integration, simply
 * steering towards static goals without avoiding any other agents or static obstacles.
 * Agents that are "selected" in the GUI will have some simple annotations that
 * show how the spatial database and engine interface can be used.
 *
 * This class is instantiated when the engine calls SimpleAIModule::createAgent().
 *
 */
class SimpleAgent : public SteerLib::AgentInterface
{
public:
	SimpleAgent();
	~SimpleAgent();
	void reset(const SteerLib::AgentInitialConditions & initialConditions, SteerLib::EngineInterface * engineInfo);
	void updateAI(float timeStamp, float dt, unsigned int frameNumber);
	void draw();

	bool enabled() { return _enabled; }
	Util::Point position() { return _position; }
	Util::Vector forward() { return _forward; }
	float radius() { return _radius; }
	const SteerLib::AgentGoalInfo & currentGoal() { return _goalQueue.front(); }
	void addGoal(const SteerLib::AgentGoalInfo & newGoal) { throw Util::GenericException("addGoals() not implemented yet for SimpleAgent"); }
	void clearGoals() { throw Util::GenericException("clearGoals() not implemented yet for SimpleAgent"); }

	/// @name The SteerLib::SpatialDatabaseItemInterface
	/// @brief These functions are required so that the agent can be used by the SteerLib::GridDatabase2D spatial database;
	/// The Util namespace helper functions do the job nicely for basic circular agents.
	//@{
	bool intersects(const Util::Ray &r, float &t) { return Util::rayIntersectsCircle2D(_position, _radius, r, t); }
	bool overlaps(const Util::Point & p, float radius) { return Util::circleOverlapsCircle2D( _position, _radius, p, radius); }
	float computePenetration(const Util::Point & p, float radius) { return Util::computeCircleCirclePenetration2D( _position, _radius, p, radius); }
	//@}


protected:
	/// Updates position, velocity, and orientation of the agent, given the force and dt time step.
	void _doEulerStep(const Util::Vector & steeringDecisionForce, float dt);

	bool _enabled;
	Util::Point _position;
	Util::Vector _velocity;
	Util::Vector _forward; // normalized version of velocity
	float _radius;
	std::queue<SteerLib::AgentGoalInfo> _goalQueue;
};

#endif
