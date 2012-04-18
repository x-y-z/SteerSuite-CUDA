//
// Copyright (c) 2009-2010 Shawn Singh, Mubbasir Kapadia, Petros Faloutsos, Glenn Reinman
// See license.txt for complete license.
//

/// @file GridDatabase2D.cpp
/// @brief Implements the SteerLib::GridDatabase2D spatial database.

#include <set>
#include <iostream>
#include <algorithm>

#include "util/GenericException.h"
#include "util/Geometry.h"
#include "util/DrawLib.h"
#include "util/Color.h"
#include "util/Misc.h"
#include "mersenne/MersenneTwister.h"

#include "interfaces/AgentInterface.h"
#include "griddatabase/GridDatabase2D.h"
#include "griddatabase/GridDatabasePlanningDomain.h"


using namespace std;
using namespace SteerLib;
using namespace Util;


//
// constructor for grid database - takes the bounds and desired number of cells
//
GridDatabase2D::GridDatabase2D(float xmin, float xmax, float zmin, float zmax, unsigned int numXCells, unsigned int numZCells, unsigned int maxItemsPerCell)
{
	if (xmin > xmax) swap(xmin, xmax);
	if (zmin > zmax) swap(xmin, xmax);

	_xOrigin = xmin;
	_zOrigin = zmin;
	_xGridSize = xmax - xmin;
	_zGridSize = zmax - zmin;
	_xInvGridSize = 1.0f / _xGridSize;
	_zInvGridSize = 1.0f / _zGridSize;
	_xNumCells = numXCells;
	_zNumCells = numZCells;
	_xCellSize = _xGridSize / ((float)numXCells);
	_zCellSize = _zGridSize / ((float)numZCells);
	_maxItemsPerCell = maxItemsPerCell;

	_allocateDatabase();
	_planningDomain = new GridDatabasePlanningDomain(this);

	cudaItems = NULL;
	hostItems = NULL;
	cudaAgentNum = 0;
	cudaObstacleNum = 0;
}


//
// alternate constructor for grid database - specify the (minx,minz) corner, and the desired grid size and number of cells
//
GridDatabase2D::GridDatabase2D(const Point & origin2D, float xExtent, float zExtent, unsigned int numXCells, unsigned int numZCells, unsigned int maxItemsPerCell)
{
	_xOrigin = origin2D.x;
	_zOrigin = origin2D.z;
	_xGridSize = xExtent;
	_zGridSize = zExtent;
	_xInvGridSize = 1.0f / _xGridSize;
	_zInvGridSize = 1.0f / _zGridSize;
	_xNumCells = numXCells;
	_zNumCells = numZCells;
	_xCellSize = _xGridSize / ((float)numXCells);
	_zCellSize = _zGridSize / ((float)numZCells);
	_maxItemsPerCell = maxItemsPerCell;

	_allocateDatabase();
	_planningDomain = new GridDatabasePlanningDomain(this);

	cudaItems = NULL;
	hostItems = NULL;
	cudaAgentNum = 0;
	cudaObstacleNum = 0;
}


//
// destructor
//
GridDatabase2D::~GridDatabase2D()
{
	delete [] _basePtr;
	delete [] _cells;
	delete _planningDomain;

	if (cudaItems != NULL)
		cudaFree(cudaItems);

	if (hostItems != NULL)
		delete [] hostItems;
}


//
// _allocateDatabase() - protected helper function called by constructors.
//
void GridDatabase2DPrivate::_allocateDatabase()
{
	if (_maxItemsPerCell > 15) {
		cerr << "WARNING: GridDatabase2D is linearly dependent on _maxItemsPerCell." << endl;
		cerr << "         right now _maxItemsPerCell is " << _maxItemsPerCell << "," << endl;
		cerr << "         performance may degrade." << endl;
	}

	unsigned int numTotalCells = _xNumCells*_zNumCells;
	unsigned int numTotalItems = numTotalCells * _maxItemsPerCell;

	_basePtr = new SpatialDatabaseItemPtr[numTotalItems];
	_cells = new GridCell[numTotalCells];

	for (unsigned int i=0; i < numTotalCells; i++) {
		// TODO: is it OK to make the traversal cost 0.0f ?? it would be more general.  need to double-check assumptions 
		// of astar lib...  is traversal cost a fixed cost to add, or is it a multiplicative factor?
		_cells[i].init( _maxItemsPerCell, _basePtr + (i*_maxItemsPerCell), 1.0f );
	}
}


//
// _clampSpatialBoundsToIndexRange - converts a bounding box to a valid range of indices in the database.
//                                  if the bounding box is entirely outside of the space covered by the
//                                  database, returns false.
//
// because it is inline, if this function ever needs to become public, place it in the .h file instead.
//
inline bool GridDatabase2DPrivate::_clampSpatialBoundsToIndexRange(float xmin, float xmax, float zmin, float zmax, unsigned int & xMinIndex, unsigned int & xMaxIndex, unsigned int & zMinIndex, unsigned int & zMaxIndex)
{
	// clamp and convert xmin to xMinIndex
	if (xmin < _xOrigin)
		xMinIndex = 0;
	else if (xmin > _xOrigin + _xGridSize)
		return false;
	else
		xMinIndex = (unsigned int) (((xmin - _xOrigin) * _xInvGridSize) * _xNumCells);

	// clamp and convert zmin to zMinIndex
	if (zmin < _zOrigin)
		zMinIndex = 0;
	else if (zmin > _zOrigin + _zGridSize)
		return false;
	else
		zMinIndex = (unsigned int) (((zmin - _zOrigin) * _zInvGridSize) * _zNumCells);

	// clamp and convert xmax to xMaxIndex
	if (xmax >= _xOrigin + _xGridSize)
		xMaxIndex = _xNumCells-1; // subtract one because xMaxIndex is included in the spatial bounds.
	else if (xmax < _xOrigin)
		return false;
	else
		xMaxIndex = (unsigned int) (((xmax - _xOrigin) * _xInvGridSize) * _xNumCells);

	// clamp and convert zmax to zMaxIndex
	if (zmax >= _zOrigin + _zGridSize)
		zMaxIndex = _zNumCells-1;  // subtract one because zMaxIndex is included in the spatial bounds.
	else if (zmax < _zOrigin)
		return false;
	else
		zMaxIndex = (unsigned int) (((zmax  - _zOrigin) * _zInvGridSize) * _zNumCells);

	return true;
}


//
// addObject() - adds the given item to the database.  Each grid cell that overlaps
//               "newBounds" will then contain a reference to the item.
//
void GridDatabase2D::addObject( SpatialDatabaseItemPtr item, const AxisAlignedBox & newBounds )
{
	// convert the spatial bounds of the object into index bounds
	unsigned int xMinIndex, xMaxIndex, zMinIndex, zMaxIndex;
	if (_clampSpatialBoundsToIndexRange(newBounds.xmin, newBounds.xmax, newBounds.zmin, newBounds.zmax, xMinIndex, xMaxIndex, zMinIndex, zMaxIndex) == false) {
		// if we get false here, the object's bounds are completely outside the database anyway.
		return;
	}

	unsigned int cellIndex;

	// iterate over all cells that overlap the bounding box of the object
	// Note that we use a shortcut to avoid wasting computation;  instead of using getCellIndexFromGridCoords(i,j) in the inner-loop,
	// we take advantage of the fact that the grid cells are contiguous in memory, and increment cellIndex directly,
	// recomputing it only when i changes.
	for (unsigned int i=xMinIndex; i<=xMaxIndex; i++) {
		cellIndex = getCellIndexFromGridCoords(i,zMinIndex);
		for (unsigned int j=zMinIndex; j<=zMaxIndex; j++) {
			_cells[cellIndex].add(item,_maxItemsPerCell, item->getTraversalCost());
			cellIndex++;
		}
	}
}


//
// removeObject() - removes an item from the grid cells that overlap with "oldBounds"
//
void GridDatabase2D::removeObject( SpatialDatabaseItemPtr item, const AxisAlignedBox &oldBounds )
{
	// convert the spatial bounds of the object into index bounds
	unsigned int xMinIndex, xMaxIndex, zMinIndex, zMaxIndex;
	if (_clampSpatialBoundsToIndexRange(oldBounds.xmin, oldBounds.xmax, oldBounds.zmin, oldBounds.zmax, xMinIndex, xMaxIndex, zMinIndex, zMaxIndex) == false) {
		// if we get false here, the object's bounds are completely outside the database anyway.
		return;
	}

	unsigned int cellIndex;

	// iterate over all cells that overlap the bounding box of the object
	// Note that we use a shortcut to avoid wasting computation;  instead of using getCellIndexFromGridCoords(i,j) in the inner-loop,
	// we take advantage of the fact that the grid cells are contiguous in memory, and increment cellIndex directly,
	// recomputing it only when i changes.
	for (unsigned int i=xMinIndex; i<=xMaxIndex; i++) {
		cellIndex = getCellIndexFromGridCoords(i,zMinIndex);
		for (unsigned int j=zMinIndex; j<=zMaxIndex; j++) {
			_cells[cellIndex].remove(item, _maxItemsPerCell, item->getTraversalCost());
			cellIndex++;
		}
	}
}


//
// updateObject() - updates the grid cells that have a reference to the item.
//
// for now this function calls removeObject() and addObject(). but later it may be
// appropriate to optimize this even further by updating only the grid cells that
// will be changed by the update.
//
//
void GridDatabase2D::updateObject( SpatialDatabaseItemPtr item, const AxisAlignedBox & oldBounds, const AxisAlignedBox & newBounds )
{
	// TODO: make an efficient "diff" between the two bounding boxes, and only iterate over the disjoint parts.
	removeObject(item, oldBounds);
	addObject(item, newBounds);
}


//
// getItemsInRange() - the protected version uses the integer index ranges.
//
void GridDatabase2D::getItemsInRange(set<SpatialDatabaseItemPtr> & neighborList, unsigned int xMinIndex, unsigned int xMaxIndex, unsigned int zMinIndex, unsigned int zMaxIndex, SpatialDatabaseItemPtr exclude)
{
	int cellIndex;

	// iterate over all grid cells in the range,
	for (unsigned int i=xMinIndex; i<=xMaxIndex; i++) {
		cellIndex = (i * _zNumCells) + zMinIndex;
		for (unsigned int j=zMinIndex; j<=zMaxIndex; j++) {
			for (unsigned int k=0; k < _maxItemsPerCell; k++) {
				if ((_cells[cellIndex]._items[k]!=NULL) && (_cells[cellIndex]._items[k]!=exclude)) {
					neighborList.insert(_cells[cellIndex]._items[k]);
				}
			}
			cellIndex++;
		}
	}
}


//
// getItemsInRange() - simply converts the spatial bounds into index range, and then calls the private getItemsInRange().
//
void GridDatabase2D::getItemsInRange(set<SpatialDatabaseItemPtr> & neighborList, float xmin, float xmax, float zmin, float zmax, SpatialDatabaseItemPtr exclude)
{
	unsigned int xMinIndex=0, xMaxIndex=0, zMinIndex=0, zMaxIndex=0;
	_clampSpatialBoundsToIndexRange(xmin, xmax, zmin, zmax, xMinIndex, xMaxIndex, zMinIndex, zMaxIndex);
	getItemsInRange(neighborList,xMinIndex,xMaxIndex,zMinIndex,zMaxIndex,exclude);
}

//
// getItemsInVisualField()
//
void GridDatabase2D::getItemsInVisualField(set<SpatialDatabaseItemPtr> & neighborList, float xmin, float xmax, float zmin, float zmax, SpatialDatabaseItemPtr exclude, const Point & position, const Vector & facingDirection, float radiusSquared)
{
	unsigned int xMinIndex=0, xMaxIndex=0, zMinIndex=0, zMaxIndex=0;
	_clampSpatialBoundsToIndexRange(xmin, xmax, zmin, zmax, xMinIndex, xMaxIndex, zMinIndex, zMaxIndex);

	int cellIndex;
	// iterate over all grid cells in the range,
	for (unsigned int i=xMinIndex; i<=xMaxIndex; i++) {
		cellIndex = getCellIndexFromGridCoords(i,zMinIndex);
		for (unsigned int j=zMinIndex; j<=zMaxIndex; j++) {
			for (unsigned int k=0; k < _maxItemsPerCell; k++) {
				SpatialDatabaseItemPtr possiblyVisibleObject = _cells[cellIndex]._items[k];


				// ignore this object if it is actually not an object (NULL), or if we are supposed to exclude it
				if ((possiblyVisibleObject==NULL) || (possiblyVisibleObject==exclude))
					continue;

				if (possiblyVisibleObject->isAgent()) {
					// three more conditions...

					// do a search to make sure it doesnt exist in the set already, if it does exist 
					// then we don't need to consider this object any further.
					if (neighborList.find(possiblyVisibleObject) != neighborList.end()) continue;

					// (1) if the agent is outside of the radius of the visual field, then forget it
					Point hisPosition = (dynamic_cast<AgentInterface*>(possiblyVisibleObject))->position();
					Vector directionToOtherAgent = hisPosition - position;
					float distSquared = directionToOtherAgent.lengthSquared();
					if (distSquared > radiusSquared) 
						continue;

					// (2) check whether the object is actually in the cone based on our facing direction
					// TODO: we shouldn't need to normalize here, because we just want to check sign, right?
					//float cosTheta = dot(directionToOtherAgent/sqrtf(distSquared),normalize(facingDirection));
					float cosTheta = dot(directionToOtherAgent/sqrtf(distSquared),normalize(facingDirection));
					if (cosTheta < 0.0f) 
						continue;

					// (3) finally, we can do the most expensive final check - checking line-of-sight.
					// previous database did not do this here, because ray tracing routines were not possible to call in the grid DB.
					// now with the virtualized interface of database items, we can.
					if (!hasLineOfSight(position, hisPosition, possiblyVisibleObject, exclude))
						continue;
					
					// if we really got this far, that means this object really is visible, so add it to the neighborList.
					neighborList.insert(possiblyVisibleObject);
				}
				else {
					// for non-Agent items, i.e. "objects" --> obstacles and such: we assume the 
					// agent will know where such items are, even if its not directly 
					// in its visual field.  so, always add it without needing to check
					// line-of-sight or distance.
					neighborList.insert(possiblyVisibleObject);
				}

			}
			cellIndex++;
		}
	}
}

void GridDatabase2D::draw()
{
#ifdef ENABLE_GUI
	Point beginning, end; // two points of a line segment to be drawn
	Point center;
	Vector xOffset(0.0f, 0.0f, 0.0f), zOffset(0.0f, 0.0f, 0.0f);
	xOffset.x = 0.5f * _xCellSize;
	zOffset.z = 0.5f * _zCellSize;

	// draw the grid itself
	DrawLib::glColor(gGray40);
	beginning.x = _xOrigin;   beginning.y = -0.01f;   beginning.z = _zOrigin;
	end.x = _xOrigin;   end.y = -0.01f;   end.z = _zOrigin + _zGridSize;
	DrawLib::drawLine(beginning, end);
	for (unsigned int i=1; i <= _xNumCells; i++) {
		float temp = _xOrigin + (i * _xCellSize);
		beginning.x = temp;
		end.x = temp;
		DrawLib::drawLine(beginning, end);
	}
	beginning.x = _xOrigin;   beginning.z = _zOrigin;
	end.x = _xOrigin + _xGridSize;   end.z = _zOrigin;
	DrawLib::drawLine(beginning, end);
	for (unsigned int j=1; j <= _zNumCells; j++) {
		float temp = _zOrigin + (j * _zCellSize);
		beginning.z = temp;
		end.z = temp;
		DrawLib::drawLine(beginning, end);
	}

#else
	throw GenericException("GridDatabase2D::draw() cannot be called, this version of SteerLib compiled without GUI functionality.");
#endif // ifdef ENABLE_GUI
}

bool GridDatabase2D::trace(const Ray & r, float & t, SpatialDatabaseItemPtr &hitObject, SpatialDatabaseItemPtr exclude, bool excludeAgents)
{
	// 1. march through grid cells
	// 2. for each grid cell:
	//   - test all items in each grid cell
	// 3. if th intersection found lies within the current grid cell, then return true with that info
	// 4. else, go on to the next grid cell.

	// just declare all vars in the first place
	unsigned int x,z; // (x,z) are 2-d coordinates and currentBin is the 1-d index to refer to the bin being traversed
	int currentBin;
	Point center;
	float xlow, xhi, zlow, zhi, xOffset, zOffset;  // spatial info for the grid cell being traversed
	float txnear, txfar, tznear, tzfar, invRayDirx, invRayDirz;
	float maxt, mint;
	bool validIntersectionFound;

	// one-time initialization stuff here
	currentBin = getCellIndexFromLocation(r.pos.x, r.pos.z);
	if (currentBin == -1) return false;
	getGridCoordinatesFromIndex(currentBin,x,z);
	invRayDirx = 1.0f / r.dir.x;
	invRayDirz = 1.0f / r.dir.z;
	xOffset = 0.5f * _xGridSize / ((float)_xNumCells);
	zOffset = 0.5f * _zGridSize / ((float)_zNumCells);

	// clamp maxt to be within the grid
	xlow = _xOrigin;
	xhi  = _xOrigin + _xGridSize;
	zlow = _xOrigin;
	zhi  = _xOrigin + _zGridSize;
	txnear = (xlow-r.pos.x) * invRayDirx;
	txfar = (xhi -r.pos.x) * invRayDirx;
	if (txnear>txfar) swap(txnear,txfar);
	tznear = (zlow-r.pos.z) * invRayDirz;
	tzfar = (zhi -r.pos.z) * invRayDirz;
	if (tznear>tzfar) swap(tznear,tzfar);
	maxt = min(r.maxt, min(txfar, tzfar));  // txnear and tznear will be "behind" the ray, since this bounding box is the entire grid and we only reach this code if the position is inside of the grid.

	// ****** URGENT TODO: ***** really should clamp mint, too...  if you are crashing when the agent is outside the grid, check this error.
	mint = r.mint;

	// set up info for the first grid cell
	getLocationFromIndex(currentBin, center);
	xlow = center.x - xOffset;
	xhi  = center.x + xOffset;
	zlow = center.z - zOffset;
	zhi  = center.z + zOffset;
	txnear = (xlow-r.pos.x) * invRayDirx;
	txfar = (xhi -r.pos.x) * invRayDirx;
	if (txnear>txfar) swap(txnear,txfar);
	tznear = (zlow-r.pos.z) * invRayDirz;
	tzfar = (zhi -r.pos.z) * invRayDirz;
	if (tznear>tzfar) swap(tznear,tzfar);

	// march through all appropriate grid cells, in order, terminating early when an *appropriate* intersection is found
	// "appropriate" means that the t value is INSIDE of the grid cell.  otherwise there is a risk of that t value not being
	// the closest one (if another object in the next grid cell is actually closer but we dont realize it yet)
	do {
		validIntersectionFound = false;
		hitObject = NULL;
		float mostRecent_maxt = min(maxt,min(txfar,tzfar)); // this way no intersection will be valid unless it was within this grid cell

		for (unsigned int i=0; i<_maxItemsPerCell; i++) {
			if ((_cells[currentBin]._items[i] != NULL) && (_cells[currentBin]._items[i] != exclude)) {

				if ((excludeAgents) && _cells[currentBin]._items[i]->isAgent())
					continue;

				float temp_t;
				bool intersected;
				intersected = false;
				Ray tempRay;
				tempRay.initWithUnitInterval(r.pos, r.dir);
				tempRay.maxt = mostRecent_maxt;
				tempRay.mint = mint;
				intersected = _cells[currentBin]._items[i]->intersects(tempRay,temp_t);
				if ((intersected) && (temp_t < mostRecent_maxt)) {
					// found a valid intersection, set all the values appropriately
					validIntersectionFound = true;
					mostRecent_maxt = temp_t;
					t = temp_t;
					hitObject = _cells[currentBin]._items[i];
				}
			}
		}

		// if a valid intersection was found in this bin, then just return
		if (validIntersectionFound) { return true; }

		// otherwise, decide which way to move to the next grid cell
		if (txfar < tzfar) {
			// then move along the x grid cells
			if (r.dir.x < 0.0f)
				x--;
			else
				x++;
			if ((x < 0) || (x > _xNumCells)) {
				return false; // no intersection found yet, and now we're out of bounds.
			}
		}
		else {
			// then move along the z grid cells
			if (r.dir.z < 0.0f)
				z--;
			else
				z++;
			if ((z < 0) || (z > _zNumCells)) {
				return false; // no intersection found yet, and now we're out of bounds.
			}
		}

		// set up info for the next grid cell
		currentBin = getCellIndexFromGridCoords(x,z);
		getLocationFromIndex(currentBin, center);
		xlow = center.x - xOffset;
		xhi  = center.x + xOffset;
		zlow = center.z - zOffset;
		zhi  = center.z + zOffset;
		txnear = (xlow-r.pos.x) * invRayDirx;
		txfar = (xhi -r.pos.x) * invRayDirx;
		if (txnear>txfar) swap(txnear,txfar);
		tznear = (zlow-r.pos.z) * invRayDirz;
		tzfar = (zhi -r.pos.z) * invRayDirz;
		if (tznear>tzfar) swap(tznear,tzfar);

	} while ( maxt > max(txnear,tznear));  // keep doing this loop until the ray's endpoint is in the cell you just traversed

	return false;

}

bool GridDatabase2D::hasLineOfSight(const Ray & r, SpatialDatabaseItemPtr exclude1, SpatialDatabaseItemPtr exclude2)
{
	// 1. march through grid cells
	// 2. for each grid cell:
	//   - test all items in each grid cell
	// 3. if th intersection found lies within the current grid cell, then return false
	// 4. else, go on to the next grid cell.

	// just declare all vars in the first place
	unsigned int x,z; // (x,z) are 2-d coordinates and currentBin is the 1-d index to refer to the bin being traversed
	int currentBin;
	Point center;
	float xlow, xhi, zlow, zhi, xOffset, zOffset;  // spatial info for the grid cell being traversed
	float txnear, txfar, tznear, tzfar, invRayDirx, invRayDirz;
	float maxt, mint;
	bool validIntersectionFound;

	// one-time initialization stuff here
	currentBin = getCellIndexFromLocation(r.pos.x, r.pos.z);
	if (currentBin == -1) return true;
	getGridCoordinatesFromIndex(currentBin,x,z);
	invRayDirx = 1.0f / r.dir.x;
	invRayDirz = 1.0f / r.dir.z;
	xOffset = 0.5f * _xGridSize / ((float)_xNumCells);
	zOffset = 0.5f * _zGridSize / ((float)_zNumCells);

	// clamp maxt to be within the grid
	xlow = _xOrigin;
	xhi  = _xOrigin + _xGridSize;
	zlow = _xOrigin;
	zhi  = _xOrigin + _zGridSize;
	txnear = (xlow-r.pos.x) * invRayDirx;
	txfar = (xhi -r.pos.x) * invRayDirx;
	if (txnear>txfar) swap(txnear,txfar);
	tznear = (zlow-r.pos.z) * invRayDirz;
	tzfar = (zhi -r.pos.z) * invRayDirz;
	if (tznear>tzfar) swap(tznear,tzfar);
	maxt = min(r.maxt, min(txfar, tzfar));  // txnear and tznear will be "behind" the ray, since this bounding box is the entire grid and we only reach this code if the position is inside of the grid.

	// ****** URGENT TODO: ***** really should clamp mint, too...  if you are crashing when the agent is outside the grid, check this error.
	mint = r.mint;

	// set up info for the first grid cell
	getLocationFromIndex(currentBin, center);
	xlow = center.x - xOffset;
	xhi  = center.x + xOffset;
	zlow = center.z - zOffset;
	zhi  = center.z + zOffset;
	txnear = (xlow-r.pos.x) * invRayDirx;
	txfar = (xhi -r.pos.x) * invRayDirx;
	if (txnear>txfar) swap(txnear,txfar);
	tznear = (zlow-r.pos.z) * invRayDirz;
	tzfar = (zhi -r.pos.z) * invRayDirz;
	if (tznear>tzfar) swap(tznear,tzfar);

	// march through all appropriate grid cells, in order, terminating early when an *appropriate* intersection is found
	// "appropriate" means that the t value is INSIDE of the grid cell.  otherwise there is a risk of that t value not being
	// the closest one (if another object in the next grid cell is actually closer but we dont realize it yet)
	do {
		validIntersectionFound = false;
		float mostRecent_maxt = min(maxt,min(txfar,tzfar)); // this way no intersection will be valid unless it was within this grid cell

		for (unsigned int i=0; i<_maxItemsPerCell; i++) {
			if ((_cells[currentBin]._items[i] != NULL) && (_cells[currentBin]._items[i] != exclude1) 
				&& (_cells[currentBin]._items[i] != exclude2) && (_cells[currentBin]._items[i]->blocksLineOfSight())) {

				float temp_t;
				bool intersected;
				intersected = false;
				Ray tempRay;
				tempRay.initWithUnitInterval(r.pos, r.dir);
				tempRay.maxt = mostRecent_maxt;
				tempRay.mint = mint;
				intersected = _cells[currentBin]._items[i]->intersects(tempRay,temp_t);
				if ((intersected) && (temp_t < mostRecent_maxt)) {
					// found a valid intersection, set all the values appropriately
					validIntersectionFound = true;
					mostRecent_maxt = temp_t;
				}
			}
		}


		// if a valid intersection was found in this bin, then just return
		if (validIntersectionFound) { return false; }

		// otherwise, decide which way to move to the next grid cell
		if (txfar < tzfar) {
			// then move along the x grid cells
			if (r.dir.x < 0.0f)
				x--;
			else
				x++;
			if ((x < 0) || (x > _xNumCells)) {
				return true; // no intersection found yet, and now we're out of bounds.
			}
		}
		else {
			// then move along the z grid cells
			if (r.dir.z < 0.0f)
				z--;
			else
				z++;
			if ((z < 0) || (z > _zNumCells)) {
				return true; // no intersection found yet, and now we're out of bounds.
			}
		}

		// set up info for the next grid cell
		currentBin = getCellIndexFromGridCoords(x,z);
		getLocationFromIndex(currentBin, center);
		xlow = center.x - xOffset;
		xhi  = center.x + xOffset;
		zlow = center.z - zOffset;
		zhi  = center.z + zOffset;
		txnear = (xlow-r.pos.x) * invRayDirx;
		txfar = (xhi -r.pos.x) * invRayDirx;
		if (txnear>txfar) swap(txnear,txfar);
		tznear = (zlow-r.pos.z) * invRayDirz;
		tzfar = (zhi -r.pos.z) * invRayDirz;
		if (tznear>tzfar) swap(tznear,tzfar);

	} while ( maxt > max(txnear,tznear));  // keep doing this loop until the ray's endpoint is in the cell you just traversed

	return true;
}

bool GridDatabase2D::hasLineOfSight(const Point & p1, const Point & p2, SpatialDatabaseItemPtr exclude1, SpatialDatabaseItemPtr exclude2)
{
	Ray r;
	r.initWithUnitInterval(p1, p2-p1);
	return hasLineOfSight(r, exclude1, exclude2);
}


Point GridDatabase2D::randomPositionWithoutCollisions(float radius, bool excludeAgents)
{
	AxisAlignedBox aab(_xOrigin, _xOrigin + _xGridSize, 0.0f, 0.0f, _zOrigin, _zOrigin + _zGridSize);
	return randomPositionInRegionWithoutCollisions(aab, radius, excludeAgents);
}

Point GridDatabase2D::randomPositionInRegionWithoutCollisions(const AxisAlignedBox & region, float radius, bool excludeAgents)
{
	static MTRand _randomNumberGenerator(2);
	return randomPositionInRegionWithoutCollisions(region, radius, excludeAgents, _randomNumberGenerator);
}

Point GridDatabase2D::randomPositionInRegionWithoutCollisions(const AxisAlignedBox & region, float radius, bool excludeAgents,  MTRand & randomNumberGenerator)
{
	Point ret(0.0f, 0.0f, 0.0f);
	bool notFoundYet;
	unsigned int numTries = 0;
	float xspan = region.xmax - region.xmin - 2*radius;
	float zspan = region.zmax - region.zmin - 2*radius;

	do {

		ret.x = region.xmin + radius + ((float)randomNumberGenerator.rand(xspan));
		ret.y = 0.0f;
		ret.z = region.zmin + radius + ((float)randomNumberGenerator.rand(zspan));

		// assume this new point has no collisions, until we find out below
		notFoundYet = false;

		// check if ret collides with anything
		set<SpatialDatabaseItemPtr> neighbors;
		neighbors.clear();
		getItemsInRange(neighbors, ret.x - radius, ret.x + radius, ret.z - radius, ret.z + radius, NULL);

		set<SpatialDatabaseItemPtr>::iterator neighbor;
		for (neighbor = neighbors.begin(); neighbor != neighbors.end(); ++neighbor) {
			if ((excludeAgents) && ((*neighbor)->isAgent())) {
				continue;
			}
			notFoundYet = (*neighbor)->overlaps(ret, radius);
			if (notFoundYet) {
				break;
			}
		}
		numTries++;
		if (numTries > 1000) {
			cerr << "ERROR: trying too hard to find a random position in region.  The region is probably already too dense." << endl;
			if (numTries > 10000) {
				throw GenericException("Gave up trying to find a random position in region.");
			}
		}
	} while (notFoundYet);

	return ret;
}


bool GridDatabase2D::planPath(unsigned int startLocation, unsigned int goalLocation, std::stack<unsigned int> & outputPlan) { 
	BestFirstSearchPlanner<GridDatabasePlanningDomain, unsigned int> gridAStarPlanner;
	gridAStarPlanner.init(_planningDomain, INT_MAX);
	return gridAStarPlanner.computePlan(startLocation, goalLocation, outputPlan);
}

void GridDatabase2D::allocateCUDAItems(int agentNum, int obstacleNum)
{
	cudaAgentNum = agentNum;
	cudaObstacleNum = obstacleNum;

	CudaSafeCall(cudaMalloc(&cudaItems, sizeof(cuda_item)*(cudaAgentNum+cudaObstacleNum)));

	hostItems = new cuda_item[cudaAgentNum+cudaObstacleNum];
}

void GridDatabase2D::addAgentCUDA(AgentInitialConditions &agentInfo, int idx)
{
	//Util::AxisAlignedBox oldBounds(_position.x-_radius, _position.x+_radius, 0.0f, 0.0f, _position.z-_radius, _position.z+_radius);

	// initialize the agent based on the initial conditions
	hostItems[idx].type = 0;
	
	hostItems[idx]._agent._position.x = agentInfo.position.x;
	hostItems[idx]._agent._position.y = agentInfo.position.y;
	hostItems[idx]._agent._position.z = agentInfo.position.z;

	hostItems[idx]._agent._forward.x = agentInfo.direction.x;
	hostItems[idx]._agent._forward.y = agentInfo.direction.y;
	hostItems[idx]._agent._forward.z = agentInfo.direction.z;


	hostItems[idx]._agent._radius = agentInfo.radius;

	hostItems[idx]._agent._velocity.x = (agentInfo.speed * Util::normalize(agentInfo.direction)).x;
	hostItems[idx]._agent._velocity.y = (agentInfo.speed * Util::normalize(agentInfo.direction)).y;
	hostItems[idx]._agent._velocity.z = (agentInfo.speed * Util::normalize(agentInfo.direction)).z;

	// compute the "new" bounding box of the agent
	Util::AxisAlignedBox newBounds(hostItems[idx]._agent._position.x-hostItems[idx]._agent._radius, hostItems[idx]._agent._position.x+hostItems[idx]._agent._radius, 
		                           0.0f, 0.0f, hostItems[idx]._agent._position.z-hostItems[idx]._agent._radius, hostItems[idx]._agent._position.z+hostItems[idx]._agent._radius);

	hostItems[idx]._agent._newBounds = hostItems[idx]._agent._oldBounds = newBounds;
	
	hostItems[idx]._agent._goalQSize = agentInfo.goals.size();

	for (unsigned int i=0; i<agentInfo.goals.size(); i++) {
		if (agentInfo.goals[i].goalType == SteerLib::GOAL_TYPE_SEEK_STATIC_TARGET) {
			hostItems[idx]._agent._goalQueue[i].x = agentInfo.goals[i].targetLocation.x;
			hostItems[idx]._agent._goalQueue[i].y = agentInfo.goals[i].targetLocation.y;
			hostItems[idx]._agent._goalQueue[i].z = agentInfo.goals[i].targetLocation.z;

			if (agentInfo.goals[i].targetIsRandom) {
				// if the goal is random, we must randomly generate the goal.
				Util::Point aGoal = randomPositionWithoutCollisions(1.0f, true);
				agentInfo.goals[i].targetIsRandom = false;
				agentInfo.goals[i].targetLocation = aGoal;
				hostItems[idx]._agent._goalQueue[i].x = agentInfo.goals[i].targetLocation.x;
				hostItems[idx]._agent._goalQueue[i].y = agentInfo.goals[i].targetLocation.y;
				hostItems[idx]._agent._goalQueue[i].z = agentInfo.goals[i].targetLocation.z;
			}
		}
		else {
			throw Util::GenericException("Unsupported goal type; SimpleAgent only supports GOAL_TYPE_SEEK_STATIC_TARGET.");
		}
	}

	hostItems[idx]._agent._enabled = true;
	hostItems[idx]._agent._curGoal = 0;

}
void GridDatabase2D::addObstacleCUDA(const ObstacleInitialConditions &obstacleInfo, int idx)
{
	int oIdx = idx + cudaAgentNum;
	hostItems[oIdx].type = 1;

	hostItems[oIdx]._obstacle._bounds = obstacleInfo;

	
	hostItems[oIdx]._obstacle._blocksLineOfSight = (hostItems[oIdx]._obstacle._bounds.ymax > 0.7) ? true : false;
}

void GridDatabase2D::fromHostToDevice()
{
	CudaSafeCall(cudaMemcpy(cudaItems, hostItems, sizeof(cuda_item)*(cudaAgentNum + cudaObstacleNum), cudaMemcpyHostToDevice));
}
void GridDatabase2D::fromDeviceToHost()
{
	CudaSafeCall(cudaMemcpy(hostItems, cudaItems, sizeof(cuda_item)*(cudaAgentNum + cudaObstacleNum), cudaMemcpyDeviceToHost));
}

int GridDatabase2D::updateAICUDA(float currentSimulationTime, float simulatonDt, unsigned int currentFrameNumber)
{
	int allDisabledAgents = 0;
	launch_updateAICUDA(cudaItems, currentSimulationTime, simulatonDt, currentFrameNumber, cudaAgentNum, cudaObstacleNum,
		                allDisabledAgents);
	return allDisabledAgents;
}

int GridDatabase2D::updateHostAgents(std::vector<SteerLib::AgentInterface*> &agentList)
{
	cudaHostAgentInfo newInfo;

	if (agentList.size() != cudaAgentNum)
		throw GenericException("agent num does not match");

	for (int i = 0; i < agentList.size(); ++i)
	{
		if (hostItems[i].type < 0)
			continue;

		newInfo._enabled = hostItems[i]._agent._enabled;
		newInfo.newBounds.xmax = hostItems[i]._agent._newBounds.xmax;
		newInfo.newBounds.ymax = hostItems[i]._agent._newBounds.ymax;
		newInfo.newBounds.zmax = hostItems[i]._agent._newBounds.zmax;
		newInfo.newBounds.xmin = hostItems[i]._agent._newBounds.xmin;
		newInfo.newBounds.ymin = hostItems[i]._agent._newBounds.ymin;
		newInfo.newBounds.zmin = hostItems[i]._agent._newBounds.zmin;

		newInfo.oldBounds.xmax = hostItems[i]._agent._oldBounds.xmax;
		newInfo.oldBounds.ymax = hostItems[i]._agent._oldBounds.ymax;
		newInfo.oldBounds.zmax = hostItems[i]._agent._oldBounds.zmax;
		newInfo.oldBounds.xmin = hostItems[i]._agent._oldBounds.xmin;
		newInfo.oldBounds.ymin = hostItems[i]._agent._oldBounds.ymin;
		newInfo.oldBounds.zmin = hostItems[i]._agent._oldBounds.zmin;

		newInfo._forward = hostItems[i]._agent._forward;
		newInfo._position = hostItems[i]._agent._position;
		newInfo._radius = hostItems[i]._agent._radius;
		newInfo._velocity = hostItems[i]._agent._velocity;

		//may be goals
		newInfo.usedGoals = hostItems[i]._agent._usedGoal;

		agentList[i]->updateWholeAgent(newInfo);
	}

	return 0;
}