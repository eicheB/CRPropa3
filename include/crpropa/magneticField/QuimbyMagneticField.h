#ifndef CRPROPA_QUIMBYMAGNETICFIELD_H
#define CRPROPA_QUIMBYMAGNETICFIELD_H

#ifdef CRPROPA_HAVE_QUIMBY

#include "crpropa/Units.h"
#include "crpropa/magneticField/MagneticField.h"

#include "quimby/MagneticField.h"

#include "crpropa/Grid.h"

#include <stdexcept>
#include <sstream>

namespace crpropa {

/**
 @class QuimbyMagneticField
 @brief Wrapper for quimby::MagneticField
 */
class QuimbyMagneticField: public MagneticField {
	quimby::ref_ptr<quimby::MagneticField> field;
public:
	QuimbyMagneticField(quimby::ref_ptr<quimby::MagneticField> field) : field(field) {
	}
	QuimbyMagneticField(quimby::MagneticField *field) : field(field) {
	}
	Vector3d getField(const Vector3d &position) const {
		quimby::Vector3f b, r = quimby::Vector3f(position.x, position.y, position.z);
		bool isGood = field->getField(r / kpc, b);
		if (!isGood) {
			std::ostringstream str;
			str << "QuimbyMagneticField: invalid position at " << position;
			throw std::runtime_error(str.str());
		}
		return Vector3d(b.x, b.y, b.z) * gauss;
	}
};
#if 1
/**
 @class QuimbyMagneticFieldAdapter
 @brief Wrapper to use crpropa::MagneticField in Quimby
 */
class QuimbyMagneticFieldAdapter: public quimby::MagneticField {
	crpropa::ref_ptr<crpropa::MagneticField> field;
public:
	QuimbyMagneticFieldAdapter(crpropa::ref_ptr<crpropa::MagneticField> field) : field(field) {

	}

	bool getField(const quimby::Vector3f &position, quimby::Vector3f &b) const {
		crpropa::Vector3d r = crpropa::Vector3d(position.x, position.y, position.z) * crpropa::kpc;
		crpropa::Vector3d B = field->getField(r);
		b = quimby::Vector3f(B.x, B.y, B.z) / gauss;
		return true;
	}
};
#endif

class ModulatedQuimbyMagneticField: public MagneticField {
	quimby::ref_ptr<quimby::MagneticField> field;	
	ref_ptr<VectorGrid> grid;
public:
	//ModulatedQuimbyMagneticField(quimby::ref_ptr<quimby::MagneticField> field, ref_ptr<VectorGrid> grid) : field(field){ 
	//}
	//ModulatedQuimbyMagneticField(quimby::MagneticField *field, ref_ptr<VectorGrid> grid) : field(field){ 
	//}	
	ModulatedQuimbyMagneticField(quimby::ref_ptr<quimby::MagneticField> field, ref_ptr<VectorGrid> grid);
	ModulatedQuimbyMagneticField(quimby::MagneticField *field, ref_ptr<VectorGrid> grid);
	//ModulatedQuimbyMagneticField() {
	//}
	
	void setGrid(ref_ptr<VectorGrid> grid);
	void setQuimby(quimby::ref_ptr<quimby::MagneticField> field);
	ref_ptr<VectorGrid> getGrid();
	quimby::ref_ptr<quimby::MagneticField> getQuimby();
	//void setReflective(bool gridReflective, bool modGridReflective);
	double getQuimbyFieldStrength(const Vector3d &position) const;	
	Vector3d getField(const Vector3d &position) const;
};

} // namespace crpropa



#endif // CRPROPA_HAVE_QUIMBY
#endif // CRPROPA_QUIMBYMAGNETICFIELD_H
