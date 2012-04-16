//
// Copyright (c) 2009-2010 Shawn Singh, Mubbasir Kapadia, Petros Faloutsos, Glenn Reinman
// See license.txt for complete license.
//

#ifndef __STEERLIB_GRID_DATABASE_PLANNING_DOMAIN_H__
#define __STEERLIB_GRID_DATABASE_PLANNING_DOMAIN_H__

/// @file GridDatabasePlanningDomain.h
/// @brief Defines the state space interface SteerLib::GridDatabasePlanningDomain, used to plan paths in the grid database.

#include "Globals.h"
#include "griddatabase/GridDatabase2D.h"
#include "planning/BestFirstSearchPlanner.h"

namespace SteerLib {


	/**
	 * @brief The internal state space of the grid database that is provided to the BestFirstSearchPlanner.
	 *
	 * This class is implemented directly in the .h file so that most compilers can inline the functions for performance.
	 *
	 * This class should not be used directly.  Instead, use the GridDatabase2D public interface which provides
	 * path-planning functionality.
	 */
	class STEERLIB_API GridDatabasePlanningDomain {
	public:
		GridDatabasePlanningDomain(SteerLib::GridDatabase2D * spatialDatabase) : _spatialDatabase(spatialDatabase) {  }

		inline bool canBeTraversed(unsigned int index) const { return (_spatialDatabase->getTraversalCost(index) < 1000.0f); }

		inline bool isAGoalState( const unsigned int & state, const unsigned int & idealGoalState) {
			return state == idealGoalState;
		}

		inline float estimateTotalCost( const unsigned int & currentState, const unsigned int & idealGoalState, float currentg) {
			unsigned int xstart, zstart, xtarget, ztarget;
			_spatialDatabase->getGridCoordinatesFromIndex(currentState, xstart, zstart);
			_spatialDatabase->getGridCoordinatesFromIndex(idealGoalState, xtarget, ztarget);
			// NOTE diffx, and diffz are signed; without the typecasting here, can get overflow.
			int diffx = ((int)xtarget) - ((int)xstart);
			int diffz = ((int)ztarget) - ((int)zstart);
			float h = sqrtf((float)(diffx*diffx + diffz*diffz));
			return currentg + h;
		}

		inline void generateTransitions( const unsigned int & currentState, const unsigned int & previousState, const unsigned int & idealGoalState, std::vector<SteerLib::DefaultAction<unsigned int> > & transitions )
		{
			transitions.reserve(7); // there will only be 7 potential actions for any given state (the eighth one would be the previous state we came from, doesn't count)
			transitions.clear();
			unsigned int x, z;
			_spatialDatabase->getGridCoordinatesFromIndex(currentState, x, z);
			//
			// THREE conditions for each potential action:
			// 1. if x+1, x-1, z+1, or z-1 are still within proper bounds of the grid
			//     - note because these are unsigned types, we have slightly different but equivalent conditional checks.
			// 2. if it doesn't transition back to the previous state
			// 3. if it can be traversed, as told by the database itself
			// note that for diagonals, the adjacent blocks must also be traversable.
			//
			if (x+1 < _spatialDatabase->getNumCellsX()) {
				// possibly add the action that transitions to location (x+1, z)
				unsigned int n = _spatialDatabase->getCellIndexFromGridCoords(x+1,z);
				if ((n != previousState)&&(canBeTraversed(n))) {
					transitions.push_back(initAction(n,_spatialDatabase->getTraversalCost(n)));
				}

				// possibly add the action that transitions to location (x+1, z+1)
				if (z+1 < _spatialDatabase->getNumCellsZ()) {
					unsigned int n = _spatialDatabase->getCellIndexFromGridCoords(x+1,z+1);
					unsigned int nAdjacent1 = _spatialDatabase->getCellIndexFromGridCoords(x,z+1);
					unsigned int nAdjacent2 = _spatialDatabase->getCellIndexFromGridCoords(x+1,z);
					if ((n != previousState)&&(canBeTraversed(n))&&(canBeTraversed(nAdjacent1))&&(canBeTraversed(nAdjacent2))) {
						transitions.push_back(initAction(n,_spatialDatabase->getTraversalCost(n) * sqrtf(2)));
					}
				}

				// possibly add the action that transitions to location (x+1, z-1)
				if (z >= 1) {
					unsigned int n = _spatialDatabase->getCellIndexFromGridCoords(x+1,z-1);
					unsigned int nAdjacent1 = _spatialDatabase->getCellIndexFromGridCoords(x,z-1);
					unsigned int nAdjacent2 = _spatialDatabase->getCellIndexFromGridCoords(x+1,z);
					if ((n != previousState)&&(canBeTraversed(n))&&(canBeTraversed(nAdjacent1))&&(canBeTraversed(nAdjacent2))) {
						transitions.push_back(initAction(n,_spatialDatabase->getTraversalCost(n) * sqrtf(2)));
					}
				}
			}

			if (x >= 1) {
				// possibly add the action that transitions to location (x-1, z)
				unsigned int n = _spatialDatabase->getCellIndexFromGridCoords(x-1,z);
				if ((n != previousState)&&(canBeTraversed(n))) {
					transitions.push_back(initAction(n,_spatialDatabase->getTraversalCost(n)));
				}

				// possibly add the action that transitions to location (x-1, z+1)
				if (z+1 < _spatialDatabase->getNumCellsZ()) {
					unsigned int n = _spatialDatabase->getCellIndexFromGridCoords(x-1,z+1);
					unsigned int nAdjacent1 = _spatialDatabase->getCellIndexFromGridCoords(x,z+1);
					unsigned int nAdjacent2 = _spatialDatabase->getCellIndexFromGridCoords(x-1,z);
					if ((n != previousState)&&(canBeTraversed(n))&&(canBeTraversed(nAdjacent1))&&(canBeTraversed(nAdjacent2))) {
						transitions.push_back(initAction(n,_spatialDatabase->getTraversalCost(n) * sqrtf(2)));
					}
				}

				// possibly add the action that transitions to location (x-1, z-1)
				if (z >= 1) {
					unsigned int n = _spatialDatabase->getCellIndexFromGridCoords(x-1,z-1);
					unsigned int nAdjacent1 = _spatialDatabase->getCellIndexFromGridCoords(x,z-1);
					unsigned int nAdjacent2 = _spatialDatabase->getCellIndexFromGridCoords(x-1,z);
					if ((n != previousState)&&(canBeTraversed(n))&&(canBeTraversed(nAdjacent1))&&(canBeTraversed(nAdjacent2))) {
						transitions.push_back(initAction(n,_spatialDatabase->getTraversalCost(n) * sqrtf(2)));
					}
				}
			}

			// possibly add the action that transitions to location (x, z+1)
			if (z+1 < _spatialDatabase->getNumCellsZ()) {
				unsigned int n = _spatialDatabase->getCellIndexFromGridCoords(x,z+1);
				if ((n != previousState)&&(canBeTraversed(n))) {
					transitions.push_back(initAction(n,_spatialDatabase->getTraversalCost(n)));
				}
			}

			// possibly add the action that transitions to location (x, z-1)
			if (z >= 1) {
				unsigned int n = _spatialDatabase->getCellIndexFromGridCoords(x,z-1);
				if ((n != previousState)&&(canBeTraversed(n))) {
					transitions.push_back(initAction(n,_spatialDatabase->getTraversalCost(n)));
				}
			}
		}

	protected:

		inline SteerLib::DefaultAction<unsigned int> & initAction(unsigned int newState, float f) {
			_tempAction.cost = f;
			_tempAction.state = newState;
			return _tempAction;
		}


		SteerLib::GridDatabase2D * _spatialDatabase;
		SteerLib::DefaultAction<unsigned int> _tempAction;
	};



} // end namespace SteerLib;


#endif

