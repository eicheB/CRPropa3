#include "crpropa/magneticField/MagneticFieldGrid.h"

namespace crpropa {

MagneticFieldGrid::MagneticFieldGrid(ref_ptr<VectorGrid> grid) {
	setGrid(grid);
}

void MagneticFieldGrid::setGrid(ref_ptr<VectorGrid> grid) {
	this->grid = grid;
}

ref_ptr<VectorGrid> MagneticFieldGrid::getGrid() {
	return grid;
}

Vector3d MagneticFieldGrid::getField(const Vector3d &pos) const {
	return grid->interpolate(pos);
}

ModulatedMagneticFieldGrid::ModulatedMagneticFieldGrid(ref_ptr<VectorGrid> grid,
		ref_ptr<ScalarGrid> modGrid) {
	grid->setReflective(false);
	modGrid->setReflective(true);
	setGrid(grid);
	setModulationGrid(modGrid);
	setInterpol(false);
}

ModulatedMagneticFieldGrid::ModulatedMagneticFieldGrid(ref_ptr<VectorGrid> grid,
		ref_ptr<ScalarGrid> modGrid, bool i2) {
	grid->setReflective(false);
	modGrid->setReflective(true);
	setGrid(grid);
	setModulationGrid(modGrid);
	setInterpol(i2); 			// in order to use interpolate2
}

void ModulatedMagneticFieldGrid::setGrid(ref_ptr<VectorGrid> g) {
	grid = g;
}

ref_ptr<VectorGrid> ModulatedMagneticFieldGrid::getGrid() {
	return grid;
}

void ModulatedMagneticFieldGrid::setModulationGrid(ref_ptr<ScalarGrid> g) {
	modGrid = g;
}

ref_ptr<ScalarGrid> ModulatedMagneticFieldGrid::getModulationGrid() {
	return modGrid;
}

void ModulatedMagneticFieldGrid::setReflective(bool gridReflective,
		bool modGridReflective) {
	grid->setReflective(gridReflective);
	modGrid->setReflective(modGridReflective);
}

// in order to be able to use interpolate2:
void ModulatedMagneticFieldGrid::setInterpol(bool i2){
	interpol2 = i2;
}
bool ModulatedMagneticFieldGrid::getInterpol() const {
	return interpol2;
}

Vector3d ModulatedMagneticFieldGrid::getField(const Vector3d &pos) const {
	/*float m = modGrid->interpolate(pos);
	Vector3d b = grid->interpolate(pos);
	return b * m;*/
	float m;	
	if (interpol2 == true) {
	m = modGrid->interpolate2(pos);
	}
	else {
	m = modGrid->interpolate(pos);	
	}
	Vector3d b = grid->interpolate(pos);
	return b * m;
}

} // namespace crpropa
