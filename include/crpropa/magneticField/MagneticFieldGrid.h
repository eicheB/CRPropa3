#ifndef CRPROPA_MAGNETICFIELDGRID_H
#define CRPROPA_MAGNETICFIELDGRID_H

#include "crpropa/magneticField/MagneticField.h"
#include "crpropa/Grid.h"

namespace crpropa {

/**
 @class MagneticFieldGrid
 @brief Magnetic field on a periodic (or reflective), cartesian grid with trilinear interpolation.

 This class wraps a VectorGrid to serve as a MagneticField.
 */
class MagneticFieldGrid: public MagneticField {
	ref_ptr<VectorGrid> grid;
public:
	MagneticFieldGrid(ref_ptr<VectorGrid> grid);
	void setGrid(ref_ptr<VectorGrid> grid);
	ref_ptr<VectorGrid> getGrid();
	Vector3d getField(const Vector3d &position) const;
};

/**
 @class MagneticFieldGrid
 @brief Modulated magnetic field on a periodic grid.

 This class wraps a VectorGrid to serve as a MagneticField.
 The field is modulated on-the-fly with a ScalarGrid.
 The VectorGrid and ScalarGrid do not need to share the same origin, spacing or size.
 */
class ModulatedMagneticFieldGrid: public MagneticField {
	ref_ptr<VectorGrid> grid;
	ref_ptr<ScalarGrid> modGrid;
	bool interpol2;						// in order to use interpolate2
public:
	ModulatedMagneticFieldGrid() {
	}
	ModulatedMagneticFieldGrid(ref_ptr<VectorGrid> grid, ref_ptr<ScalarGrid> modGrid);
	ModulatedMagneticFieldGrid(ref_ptr<VectorGrid> grid, ref_ptr<ScalarGrid> modGrid, bool interpol2); // in order to use interpolate2
	void setGrid(ref_ptr<VectorGrid> grid);
	void setModulationGrid(ref_ptr<ScalarGrid> modGrid);
	ref_ptr<VectorGrid> getGrid();
	ref_ptr<ScalarGrid> getModulationGrid();
	void setReflective(bool gridReflective, bool modGridReflective);
	Vector3d getField(const Vector3d &position) const;
	void setInterpol(bool interpol2);			// in order to use interpolate2 
	bool getInterpol() const; 				// in order to use interpolate2
};

} // namespace crpropa

#endif // CRPROPA_MAGNETICFIELDGRID_H
