//
// Copyright (c) 2009-2010 Shawn Singh, Mubbasir Kapadia, Petros Faloutsos, Glenn Reinman
// See license.txt for complete license.
//

/// @file CommandLineEngineDriver.cpp
/// @brief Implements the CommandLineEngineDriver functionality.
///
/// @todo
///   - update documentation in this file
///

#include <iostream>
#include "SteerLib.h"
#include "core/CommandLineEngineDriver.h"

using namespace std;
using namespace SteerLib;
using namespace Util;

//
// constructor
//
CommandLineEngineDriver::CommandLineEngineDriver()
{
	_alreadyInitialized = false;
	_engine = NULL;
}


//
// init()
//
void CommandLineEngineDriver::init(SteerLib::SimulationOptions * options)
{
	if (_alreadyInitialized) {
		throw GenericException("CommandLineEngineDriver::init() - should not call this function twice.\n");
	}

	_alreadyInitialized = true;

	_engine = new SimulationEngine();
	_engine->init(options, this);
}


//
// run() - never returns, will end the program properly when appropriate
//
void CommandLineEngineDriver::run()
{
	bool verbose = true;  // TODO: make this a user option...
	bool done = false;

	if (verbose) std::cout << "\rInitializing...\n";
	_engine->initializeSimulation();

	if (verbose) std::cout << "\rPreprocessing...\n";
	_engine->preprocessSimulation();

	// loop until the engine tells us its done
	while (!done) {
		if (verbose) std::cout << "\rFrame Number:   " << _engine->getClock().getCurrentFrameNumber();
		done = !_engine->update(false);
	}

	if (verbose) std::cout << "\rFrame Number:   " << _engine->getClock().getCurrentFrameNumber() << std::endl;

	if (verbose) std::cout << "\rPostprocessing...\n";
	_engine->postprocessSimulation();

	if (verbose) std::cout << "\rCleaning up...\n";
	_engine->cleanupSimulation();

	if (verbose) std::cout << "\rDone.\n";

}

//
// finish() - cleans up.
//
void CommandLineEngineDriver::finish()
{
	_engine->finish();

	delete _engine;
}

