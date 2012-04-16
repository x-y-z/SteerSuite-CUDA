//
// Copyright (c) 2009-2010 Shawn Singh, Mubbasir Kapadia, Petros Faloutsos, Glenn Reinman
// See license.txt for complete license.
//

/// @file SimpleAIModule.cpp
/// @brief Implements the SimpleAIModule plugin.


#include "SteerLib.h"
#include "SimulationPlugin.h"
#include "SimpleAIModule.h"
#include "SimpleAgent.h"


// globally accessible to the simpleAI plugin
SteerLib::EngineInterface * gEngine;
SteerLib::GridDatabase2D * gSpatialDatabase;



PLUGIN_API SteerLib::ModuleInterface * createModule()
{
	return new SimpleAIModule;
}

PLUGIN_API void destroyModule( SteerLib::ModuleInterface*  module )
{
	delete module;
}


void SimpleAIModule::init( const SteerLib::OptionDictionary & options, SteerLib::EngineInterface * engineInfo )
{
	gEngine = engineInfo;
	gSpatialDatabase = engineInfo->getSpatialDatabase();
}

void SimpleAIModule::finish()
{
	// nothing to do here
}

SteerLib::AgentInterface * SimpleAIModule::createAgent()
{
	return new SimpleAgent; 
}

void SimpleAIModule::destroyAgent( SteerLib::AgentInterface * agent )
{
	delete agent;
}
