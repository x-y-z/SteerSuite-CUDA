//
// Copyright (c) 2009-2010 Shawn Singh, Mubbasir Kapadia, Petros Faloutsos, Glenn Reinman
// See license.txt for complete license.
//

/// @file TestCasePlayerModule.cpp
/// @brief Implements the TestCasePlayerModule built-in module.

#include "modules/TestCasePlayerModule.h"
#include "testcaseio/TestCaseIO.h"
#include "util/Misc.h"

using namespace SteerLib;

void TestCasePlayerModule::init( const SteerLib::OptionDictionary & options, SteerLib::EngineInterface * engineInfo ) {
	_engine = engineInfo;
	_testCaseFilename = "";
	_aiModuleName = "";
	_aiModuleSearchPath = "";
	_aiModule = NULL;
	_obstacles.clear();

	// parse command line options
	SteerLib::OptionDictionary::const_iterator optionIter;
	for (optionIter = options.begin(); optionIter != options.end(); ++optionIter) {
		if ((*optionIter).first == "testcase") {
			_testCaseFilename = (*optionIter).second;
		}
		else if ((*optionIter).first == "ai") {
			_aiModuleName = (*optionIter).second;
		}
		else {
			throw Util::GenericException("unrecognized option \"" + Util::toString((*optionIter).first) + "\" given to testCasePlayer module.");
		}
	}

	// load the module if it is not already
	if ( !_engine->isModuleLoaded(_aiModuleName) ) {
		_engine->loadModule(_aiModuleName,_aiModuleSearchPath,"");
		/// @todo add a boolean that flags this module was loaded here, so that we can unload the module after the simulation is done.
	}

	// get a pointer to the loaded module
	_aiModule = _engine->getModule(_aiModuleName);


#ifdef ENABLE_GUI
#ifdef ENABLE_QT

	QMainWindow * mainWindow = (QMainWindow*)(_engine->getEngineController()->getQtMainWindow());

	if(mainWindow != NULL) {

		/// @todo why is this tooltip attribute here, and not in the engine driver?
		mainWindow->QWidget::setAttribute(Qt::WA_AlwaysShowToolTips);

		_testCasePlayerDockWidget = new QDockWidget("Test Case Player");
		_testCasePlayerWidget = new SteerSimQt::TestCasePlayerWidget(this, _engine);
		_testCasePlayerDockWidget->setWidget(_testCasePlayerWidget);
		mainWindow->addDockWidget(Qt::RightDockWidgetArea, _testCasePlayerDockWidget);
	}
	else {
		_testCasePlayerDockWidget = NULL;
		_testCasePlayerWidget = NULL;
	}
#endif
#endif
}

void TestCasePlayerModule::initializeSimulation() {

	SteerLib::TestCaseReader * testCaseReader;

	std::string testCasePath;

	// try to find the test case in several ways:
	if (Util::fileCanBeOpened(_testCaseFilename)) {
		testCasePath = _testCaseFilename;
	}
	else if (Util::fileCanBeOpened( _engine->getTestCaseSearchPath() + _testCaseFilename )) {
		testCasePath = _engine->getTestCaseSearchPath() + _testCaseFilename;
	}
	else if (Util::fileCanBeOpened(_testCaseFilename + ".xml")) {
		testCasePath = _testCaseFilename + ".xml";
	}
	else if (Util::fileCanBeOpened( _engine->getTestCaseSearchPath() + _testCaseFilename + ".xml" )) {
		testCasePath = _engine->getTestCaseSearchPath() + _testCaseFilename + ".xml";
	}
	else {
		throw Util::GenericException("Could not find test case " + _testCaseFilename + ".");
	}

	// open the test case
	testCaseReader = new SteerLib::TestCaseReader();
	testCaseReader->readTestCaseFromFile(testCasePath);


	for (unsigned int i=0; i < testCaseReader->getNumObstacles(); i++) {
		const SteerLib::ObstacleInitialConditions & ic = testCaseReader->getObstacleInitialConditions(i);
		SteerLib::BoxObstacle * b;
		b = new SteerLib::BoxObstacle(ic.xmin, ic.xmax, ic.ymin, ic.ymax, ic.zmin, ic.zmax);
		_obstacles.push_back(b);
		_engine->addObstacle(b);
		_engine->getSpatialDatabase()->addObject( b, b->getBounds());
	}


	for (unsigned int i=0; i < testCaseReader->getNumAgents(); i++) {
		const SteerLib::AgentInitialConditions & ic = testCaseReader->getAgentInitialConditions(i);
		_engine->createAgent( ic, _aiModule );
	}


	delete testCaseReader;


#ifdef ENABLE_GUI
#ifdef ENABLE_QT
	if (_testCasePlayerWidget != NULL) _testCasePlayerWidget->setFrameNumber(_engine->getClock().getCurrentFrameNumber());
#endif
#endif

}

void TestCasePlayerModule::postprocessFrame(float timeStamp, float dt, unsigned int frameNumber) {
#ifdef ENABLE_GUI
#ifdef ENABLE_QT
	if (_testCasePlayerWidget != NULL) _testCasePlayerWidget->setFrameNumber(frameNumber);
#endif
#endif
}


void TestCasePlayerModule::cleanupSimulation() {
	for (unsigned int i=0; i<_obstacles.size(); i++) {
		_engine->getSpatialDatabase()->removeObject( _obstacles[i], _obstacles[i]->getBounds());
		_engine->removeObstacle(_obstacles[i]);
		delete _obstacles[i];
	}
	_obstacles.clear();

	_engine->destroyAllAgentsFromModule(_aiModule);
}


void TestCasePlayerModule::finish() {

#ifdef ENABLE_GUI
#ifdef ENABLE_QT

	QMainWindow * mainWindow = (QMainWindow*)(_engine->getEngineController()->getQtMainWindow());
	if(mainWindow != NULL) {
		mainWindow->removeDockWidget(_testCasePlayerDockWidget);
		delete _testCasePlayerDockWidget;
		delete _testCasePlayerWidget; // may not be necessary, do Qt automatically take control and de-allocate children widgets?
	}

#endif
#endif

}
