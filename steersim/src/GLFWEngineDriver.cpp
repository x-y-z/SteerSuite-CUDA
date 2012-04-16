//
// Copyright (c) 2009-2010 Shawn Singh, Mubbasir Kapadia, Petros Faloutsos, Glenn Reinman
// See license.txt for complete license.
//

/// @file GLFWEngineDriver.cpp
/// @brief Implements the GLFWEngineDriver functionality.
///
/// @todo
///   - update documentation in this file
///

#ifdef ENABLE_GUI
#ifdef ENABLE_GLFW

#include <iostream>
#include "SteerLib.h"
#include "core/GLFWEngineDriver.h"

#include "glfw/include/GL/glfw.h"

using namespace std;
using namespace SteerLib;
using namespace Util;

#define MULTISAMPLE_ARB 0x809D

//
// callback wrappers for GLFW
//
static void GLFWCALL processWindowResizedEvent(int width, int height) { GLFWEngineDriver::getInstance()->processWindowResizedEvent(width,height); }
static void GLFWCALL processKeyPressEvent(int key, int action) { GLFWEngineDriver::getInstance()->processKeyPressEvent(key, action); }
static void GLFWCALL processMouseButtonEvent(int button, int action) { GLFWEngineDriver::getInstance()->processMouseButtonEvent(button, action); }
static void GLFWCALL processMouseMovementEvent(int x, int y) { GLFWEngineDriver::getInstance()->processMouseMovementEvent(x,y); }
static void GLFWCALL processMouseWheelEvent(int pos) { GLFWEngineDriver::getInstance()->processMouseWheelEvent(pos); }
static int GLFWCALL processWindowCloseEvent() { return GLFWEngineDriver::getInstance()->processWindowCloseEvent(); }


//
// getInstance()
//
// Singleton trick:  with the static instance in this function, we are guaranteed that the 
// constructor will be called before the instance is ever retrieved.
//
GLFWEngineDriver * GLFWEngineDriver::getInstance()
{
	static GLFWEngineDriver * singletonInstance = new GLFWEngineDriver();
	return singletonInstance;
}


//
// constructor
//
GLFWEngineDriver::GLFWEngineDriver()
{
	// Note that many of these values will be changed during init() anyway
	_alreadyInitialized = false;
	_engine = NULL;
	_mouseX = 0;
	_mouseY = 0;
	_wheelPos = 0;
	_moveCameraOnMouseMotion = false;
	_rotateCameraOnMouseMotion = false;
	_zoomCameraOnMouseMotion = false;
	_multisampleAntialiasingSupported = false; // until we find out later
	_agentNearestToMouse = NULL;
	_nextScreenshotNumber = 0;
	_dumpFrames = false;
	_done = false;

	_useAntialiasing = true;
	_canUseMouseToSelectAgents = true;
	_canUseMouseWheelZoom = true;

	_options = NULL;
}


//
// init()
//
void GLFWEngineDriver::init(SimulationOptions * options)
{
	if (_alreadyInitialized) {
		throw GenericException("GLFWEngineDriver::init() was called twice, but it should only be called once.");
	}

	_options = options;

	_alreadyInitialized = true;
	_done = false;
	_dumpFrames = false;

	_engine = new SimulationEngine();
	_engine->init(_options, this);


	_useAntialiasing = _options->guiOptions.useAntialiasing;
	_canUseMouseToSelectAgents = _options->guiOptions.canUseMouseSelection;
	_canUseMouseWheelZoom = _options->guiOptions.canUseMouseWheelZoom;
	_paused = _options->glfwEngineDriverOptions.pausedOnStart;

	_initGLFW();
	_checkGLCapabilities();
	_initGL(); // calls _engine->initGL() among other things...

	DrawLib::init();
}

void GLFWEngineDriver::finish()
{
	_engine->finish();

	delete _engine;

	glfwTerminate();
}


void GLFWEngineDriver::_initGLFW()
{
	// initialize glfw
	if (glfwInit() == GL_FALSE) {
		throw GenericException("Initializing GLFW failed.\n");
	}

	// specify some configuration that needs to happen before creating the window
	glfwOpenWindowHint( GLFW_FSAA_SAMPLES, 4 );

	// create the glfw window
	if (glfwOpenWindow( _options->glfwEngineDriverOptions.windowSizeX, _options->glfwEngineDriverOptions.windowSizeY, 8, 8, 8, 8, 24, 8, GLFW_WINDOW) == GL_FALSE) {
		throw GenericException("Could not open a window using glfwOpenWindow().");
	}

	// specify some configuration that needs to happen after creating the window
	glfwSetWindowTitle( _options->glfwEngineDriverOptions.windowTitle.c_str() );
	glfwSetWindowPos( _options->glfwEngineDriverOptions.windowPositionX, _options->glfwEngineDriverOptions.windowPositionY);
	int interval = _options->guiOptions.useVsync ? 1 : 0;
	glfwSwapInterval(interval);

	// register GLFW callbacks
	// note the "::" is needed to resolve the static global functions, not the GLFWEngineDriver member functions.
	glfwSetWindowSizeCallback( ::processWindowResizedEvent );
	glfwSetKeyCallback( ::processKeyPressEvent );
	glfwSetMouseButtonCallback( ::processMouseButtonEvent );
	glfwSetMousePosCallback( ::processMouseMovementEvent );
	glfwSetMouseWheelCallback( ::processMouseWheelEvent );
	glfwSetWindowCloseCallback( ::processWindowCloseEvent );

}

/// @todo 
///   - need to properly initialize the antialiasing mode based on user options
///     and initialize _useAntialising at the same time.
///   - figure out how to approrpiately toggle vsync
void GLFWEngineDriver::_checkGLCapabilities()
{
	if (glfwExtensionSupported("GL_ARB_multisample") == GL_TRUE) {
		_multisampleAntialiasingSupported = true;
	}
	else {
		_multisampleAntialiasingSupported = false;
	}


}

//
// initGL()
//
void GLFWEngineDriver::_initGL()
{
	_engine->initGL();

	if (_multisampleAntialiasingSupported && _useAntialiasing) {
		glEnable( MULTISAMPLE_ARB );
	} else {
		glDisable( MULTISAMPLE_ARB );
	}

	int w,h;
	glfwGetWindowSize( &w, &h );
	processWindowResizedEvent(w,h);
}


//
// run() - returns if the simulation is done, but there are other ways the program will exit (e.g. user closes the window)
//
void GLFWEngineDriver::run()
{
	bool verbose = true;  // TODO: make this a user option ??.

	if (verbose) std::cout << "\rInitializing...\n";
	_engine->initializeSimulation();

	if (verbose) std::cout << "\rPreprocessing...\n";
	_engine->preprocessSimulation();

	if (verbose) std::cout << "\rSimulation is running...\n";
	while (!_done) {

		// Finding the agent closest to mouse uses an exhaustive search across all agents.
		// the algorithm could probably be better (i.e. use the spatial database to find it)
		// but its not (yet) worth implementing that way.  Just change DEFAULT_CAN_USE_MOUSE_SELECTION to 
		// false if you want to disable it.
		if (_canUseMouseToSelectAgents) {
			_findClosestAgentToMouse();
		}

		// Update the AI.
		if (_engine->update(_paused) == false) {
			// The engine indicated the simulaton should finish
			_done = true;
		}
		else {
			// The simulation is continuing, so draw everything using openGL.
			_drawScene();
		}
	}

	if (verbose) std::cout << "\rPostprocessing...\n";
	_engine->postprocessSimulation();

	if (verbose) std::cout << "\rCleaning up...\n";
	_engine->cleanupSimulation();

	if (verbose) std::cout << "\rDone.\n";
}

void GLFWEngineDriver::_drawScene()
{
	// get the camera from the engine
	Camera & cam = _engine->getCamera();

	// clear color and depth buffers
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	cam.apply();
	DrawLib::positionLights();
	_engine->draw();


	// check for errors
	int error = glGetError();
	if (error != GL_NO_ERROR) {
		std::cerr << "An OpenGL error occurred: " << gluErrorString(error) << "\n";
		throw GenericException("OpenGL error occurred in GLFWEngineDriver::_drawScene().");
	}

	if (_dumpFrames) writeScreenCapture();

	// double buffering, swap back and front buffers
	//glFlush();
	glfwSwapBuffers();
}


//
// drawGUI() - draw any GUI elements or global GUI annotations that have nothing to do with the simulation.
//
void GLFWEngineDriver::_drawGUI()
{
	if (_agentNearestToMouse != NULL) DrawLib::drawHighlight(_agentNearestToMouse->position(), _agentNearestToMouse->forward(), _agentNearestToMouse->radius());
}



void GLFWEngineDriver::writeScreenCapture()
{
	cerr << "WARNING: no screenshot taken, feature not implemented for GLFW yet.\n";
}



void GLFWEngineDriver::_findClosestAgentToMouse()
{
	// gotta find mouse/screen position in world coords

	// Since screen is 2d and world coords is 3d, we convert screen
	// position into a ray originating from the camera position and
	// going in the direction indicated by the mouse pos

	// next bit is apparently standard openGl procedure to getting what
	// we want:

	double modelView[16];
	double projection[16];
	int viewport[4];

	double x1, y1, z1, x2, y2, z2;

	glGetDoublev(GL_MODELVIEW_MATRIX, modelView);
	glGetDoublev(GL_PROJECTION_MATRIX, projection);
	glGetIntegerv(GL_VIEWPORT, viewport);

	int x = _mouseX;
	int y = viewport[3] - _mouseY;
	gluUnProject(x, y, 0, modelView, projection, viewport, &x1, &y1, &z1);
	gluUnProject(x, y, 1, modelView, projection, viewport, &x2, &y2, &z2);

	// Our ray: originates at the camara and goes somewhere
	Vector cameraRay = Vector((float)(x2-x1), (float)(y2-y1), (float)(z2-z1));
	Point cameraPos = _engine->getCamera().position();

	// to determine which agents are closest to the mouse click, figure out which
	// agents are closest to the ray
	// NOTE: assuming that agents are all at the y=0 plane.
	float t = -cameraPos.y / cameraRay.y;  // this gives us the t value for where the point along the ray would be on the y=0 plane.
	Point locationOfMouseOnYPlane = cameraPos + t*cameraRay;

	// nearest: nearest agent to mouse pos
	// nearestAndContained: nearest agent to mouse pos and the point mouse is within radius

	AgentInterface* nearest = NULL;

	float distanceToNearest = FLT_MAX;

	const std::vector<AgentInterface*> & agents = _engine->getAgents();
	for(std::vector<AgentInterface*>::const_iterator i = agents.begin(); i != agents.end(); ++i)
	{
		float dist = distanceBetween((*i)->position(), locationOfMouseOnYPlane);
		// if its the closest one, but also within some distance threshold
		if ((dist < (*i)->radius()+0.5f) && (dist < distanceToNearest)) {
			nearest = (*i);
			distanceToNearest = dist;
		}
	}
	_agentNearestToMouse = nearest;

}


//=====================================================
// GLUT callback functions
//=====================================================


void GLFWEngineDriver::processWindowResizedEvent(int width, int height)
{
	_engine->resizeGL(width, height);
}


void GLFWEngineDriver::processKeyPressEvent(int key, int action)
{

	if ((key == _options->keyboardBindings.quit) && (action==GLFW_PRESS)) {
		_done = true;
	}
	else if ((key == _options->keyboardBindings.toggleAntialiasing) && (action==GLFW_PRESS)) {
		if (!_multisampleAntialiasingSupported) {
			cerr << "WARNING: toggling antialiasing may have no effect; antialiasing does not seem to be supported.\n";
		}
		_useAntialiasing = !_useAntialiasing;
		if (_useAntialiasing) {
			glEnable( MULTISAMPLE_ARB );
		} else {
			glDisable( MULTISAMPLE_ARB );
		}
		std::cout << "Antialiasing turned " << (_useAntialiasing ? "on" : "off") << std::endl;
	}
	else if ((key == _options->keyboardBindings.printCameraInfo) && (action==GLFW_PRESS)) {
		cout << "CAMERA INFO:" << endl;
		cout << "  Position:     " << _engine->getCamera().position() << endl;
		cout << "  LookAt:       " << _engine->getCamera().lookat() << endl;
		cout << "  Up:           " << _engine->getCamera().up() << endl;
	} 
	else if ((key == _options->keyboardBindings.stepForward) && (action==GLFW_PRESS)) {
		pauseAndStepOneFrame();
	}
	else if ((key == _options->keyboardBindings.pause) && (action==GLFW_PRESS)) {
		togglePausedState();
	}
	else if ((key == (int)'F') && (action==GLFW_PRESS)) {
		cout << "Frame Number " << _engine->getClock().getCurrentFrameNumber() << "\n";
	}
	else if ((key == _options->keyboardBindings.takeScreenshot) && (action==GLFW_PRESS)) {
		writeScreenCapture();
	}
	else if ((key == _options->keyboardBindings.startDumpingFrames) && (action==GLFW_PRESS)) {
		_dumpFrames = true;
	}
	else if ((key == _options->keyboardBindings.stopDumpingFrames) && (action==GLFW_PRESS)) {
		_dumpFrames = false;
	}
	else {
		if (_engine) _engine->processKeyboardInput( key, action);
	}
}


void GLFWEngineDriver::processMouseButtonEvent(int button, int action)
{
	bool controlKeyPressed = (  (glfwGetKey(GLFW_KEY_LCTRL)==GLFW_PRESS) || (glfwGetKey(GLFW_KEY_RCTRL)==GLFW_PRESS) );

	if ((button ==  _options->mouseBindings.selectAgent) && (action == GLFW_PRESS)) {
		if(!controlKeyPressed) {
			if (_agentNearestToMouse != NULL) {
				if (!_engine->isAgentSelected(_agentNearestToMouse)) {
					_engine->selectAgent(_agentNearestToMouse);
					unsigned int i;
					for (i=0; i< _engine->getAgents().size(); i++) {
						if ( _engine->getAgents()[i] == _agentNearestToMouse ) {
							break;
						}
					}
					if (_agentNearestToMouse != NULL) cerr << "selected agent #" << i << " (total " << _engine->getSelectedAgents().size() << " agents are currently selected)." << endl;
				}
				else {
					_engine->unselectAgent(_agentNearestToMouse);
					unsigned int i;
					for (i=0; i< _engine->getAgents().size(); i++) {
						if ( _engine->getAgents()[i] == _agentNearestToMouse ) {
							break;
						}
					}
					if (_agentNearestToMouse != NULL) cerr << "un-selected agent #" << i << " (total " << _engine->getSelectedAgents().size() << " agents are currently selected)." << endl;
				}
			}
			else {
				_engine->unselectAllAgents();
			}
		}
	}

	if ((button == _options->mouseBindings.moveCamera) && (action == GLFW_PRESS)) {
		if(controlKeyPressed) _moveCameraOnMouseMotion = true;
	}

	if ((button == _options->mouseBindings.moveCamera) && (action == GLFW_RELEASE)) {
		_moveCameraOnMouseMotion = false;
	}

	if ((button == _options->mouseBindings.rotateCamera) && (action == GLFW_PRESS)) {
		if(controlKeyPressed) _rotateCameraOnMouseMotion = true;
	}

	if ((button == _options->mouseBindings.rotateCamera) && (action == GLFW_RELEASE)) {
		_rotateCameraOnMouseMotion = false;
	}

	if ((button == _options->mouseBindings.zoomCamera) && (action == GLFW_PRESS)) {
		if(controlKeyPressed) _zoomCameraOnMouseMotion = true;
	}

	if ((button == _options->mouseBindings.zoomCamera) && (action == GLFW_RELEASE)) {
		_zoomCameraOnMouseMotion = false;
	}
}


void GLFWEngineDriver::processMouseMovementEvent(int x, int y)
{

	// get mouse changes
	int deltaX = x - _mouseX;
	int deltaY = y - _mouseY;

	// update mouse position
	_mouseX = x;
	_mouseY = y;

	// camera rotate
	if(_rotateCameraOnMouseMotion)
	{
		float xAdjust = -deltaX * _options->guiOptions.mouseRotationFactor;
		float yAdjust = deltaY * _options->guiOptions.mouseRotationFactor;

		_engine->getCamera().nudgeRotate(yAdjust, xAdjust);
	}

	// camera zoom
	if(_zoomCameraOnMouseMotion)
	{
		float yAdjust = deltaY * _options->guiOptions.mouseZoomFactor;
		_engine->getCamera().nudgeZoom(yAdjust);
	}

	// camera move
	if(_moveCameraOnMouseMotion)
	{
		float xAdjust = deltaX * _options->guiOptions.mouseMovementFactor;
		float yAdjust = deltaY * _options->guiOptions.mouseMovementFactor;

		_engine->getCamera().nudgePosition(xAdjust, yAdjust);
	}
}

void GLFWEngineDriver::processMouseWheelEvent(int pos)
{
	if (_canUseMouseWheelZoom) {
		int deltaWheel = pos - _wheelPos;
		_wheelPos = pos;
		if (deltaWheel != 0) {
			float wheelAdjust = -20.0f * deltaWheel * _options->guiOptions.mouseZoomFactor;
			_engine->getCamera().nudgeZoom(wheelAdjust);
		}
	}
}

int GLFWEngineDriver::processWindowCloseEvent()
{
	// allow the user to exit at any time
	_done = true;
	return GL_TRUE;
}

#endif // ifdef ENABLE_GLFW
#endif // ifdef ENABLE_GUI
