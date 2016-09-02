#include "crpropa/magneticField/QuimbyMagneticField.h"
#include "crpropa/magneticField/MagneticFieldGrid.h"

namespace crpropa {

ModulatedQuimbyMagneticField::ModulatedQuimbyMagneticField(quimby::ref_ptr<quimby::MagneticField> field, ref_ptr<VectorGrid> grid) : field(field) {
	setGrid(grid);
}
ModulatedQuimbyMagneticField::ModulatedQuimbyMagneticField(quimby::MagneticField *field, ref_ptr<VectorGrid> grid) : field(field) {
	setGrid(grid);
}

/*

ModulatedQuimbyMagneticField::ModulatedQuimbyMagneticField(quimby::ref_ptr<quimby::MagneticField> field,
		ref_ptr<VectorGrid> grid) {
	//grid->setReflective(false);
	//modGrid->setReflective(true);
	setGrid(grid);
	//setQuimby(field);
}
ModulatedQuimbyMagneticField::ModulatedQuimbyMagneticField(quimby::MagneticField *field,
		ref_ptr<VectorGrid> grid) {
	//grid->setReflective(false);
	//modGrid->setReflective(true);
	setGrid(grid);
	//setQuimby(field);
}
*/

void ModulatedQuimbyMagneticField::setGrid(ref_ptr<VectorGrid> g) {
	grid = g;
}

ref_ptr<VectorGrid> ModulatedQuimbyMagneticField::getGrid() {
	return grid;
}

void ModulatedQuimbyMagneticField::setQuimby(quimby::ref_ptr<quimby::MagneticField> f) {
	field = f;
}

quimby::ref_ptr<quimby::MagneticField> ModulatedQuimbyMagneticField::getQuimby() {
	return field;
}

/*
Vector3d ModulatedQuimbyMagneticField::getQuimbyFieldStrength(const Vector3d &position) const {
	quimby::Vector3f b, r = quimby::Vector3f(position.x, position.y, position.z);
	bool isGood = field->getField(r / kpc, b);
	if (!isGood) {
		std::ostringstream str;
		str << "QuimbyMagneticField: invalid position at " << position;
		throw std::runtime_error(str.str());
	}
	return Vector3d(b.x, b.y, b.z) * gauss;
}
*/
double ModulatedQuimbyMagneticField::getQuimbyFieldStrength(const Vector3d &position) const {
	quimby::Vector3f b, r = quimby::Vector3f(position.x, position.y, position.z);
	bool isGood = field->getField(r / kpc, b);
	if (!isGood) {
		std::ostringstream str;
		str << "QuimbyMagneticField: invalid position at " << position;
		throw std::runtime_error(str.str());
	}
	double Bstr = Vector3d(b.x, b.y, b.z).getR();
	//return Vector3d(b.x, b.y, b.z) * gauss;
	return Bstr * gauss;
}
/*
void ModulatedQuimbyMagneticField::setReflective(bool gridReflective,
		bool modGridReflective) {
	grid->setReflective(gridReflective);
	modGrid->setReflective(modGridReflective);
}
*/

Vector3d ModulatedQuimbyMagneticField::getField(const Vector3d &pos) const {
	double QuimbyField = getQuimbyFieldStrength(pos);
	Vector3d Turb = grid->interpolate(pos);
	return Turb * QuimbyField;
}

} // namespace crpropa
