#ifndef CRPROPA_QUIMBYMAGNETICFIELD_H
#define CRPROPA_QUIMBYMAGNETICFIELD_H

#ifdef CRPROPA_HAVE_QUIMBY

#include "crpropa/Units.h"
#include "crpropa/magneticField/MagneticField.h"

#include "quimby/MagneticField.h"

#include <stdexcept>
#include <sstream>
#include <cmath>

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
class ExtendedQuimbyMagneticField: public MagneticField {
	quimby::ref_ptr<quimby::MagneticField> field;
    const Vector3d sourcePos;
public:
	ExtendedQuimbyMagneticField(quimby::ref_ptr<quimby::MagneticField> field, const Vector3d &sourcePos) : field(field), sourcePos(sourcePos) {
	}
	ExtendedQuimbyMagneticField(quimby::MagneticField *field, const Vector3d &sourcePos) : field(field), sourcePos(sourcePos) {
	}
	Vector3d getField(const Vector3d &position) const {
		quimby::Vector3f b, r = quimby::Vector3f(position.x, position.y, position.z);
        bool isGood;
        //const Vector3d sourceDis = position - sourcePos;
        // spherical reflection:
        /*
        if (sourceDis.getR()<sourcePos.min()) {
            isGood = field->getField(r / kpc, b);
        }
		else {
            //const Vector3d positionB = 2. * sourcePos + 2.*sourcePos.min() * sourceDis / sourceDis.getR() - position;
            Vector3d sourceDisNew = sourceDis;
            Vector3d positionNew = position;
            Vector3d positionB;
            while (sourceDisNew.getR()>=sourcePos.min()){
                positionB = 2. * sourcePos + 2.*sourcePos.min() * sourceDisNew / sourceDisNew.getR() - positionNew;
                positionNew = positionB;
                sourceDisNew = positionB - sourcePos;
            }
            const Vector3d NewsourceDis = positionB - sourcePos;
            if (NewsourceDis.getR()>sourcePos.min()) {
                isGood = field->getField(r / kpc, b);
            }
            r = quimby::Vector3f(positionB.x, positionB.y, positionB.z);
            isGood = field->getField(r / kpc, b);
		}
        */
        // cubic reflection:
        //if (sourceDis.getR()<sourcePos.min()) {
        double cos45 = 1./sqrt(2.);
        if (position.x - sourcePos.x * (1.-cos45) > 0. && sourcePos.x * (1.+cos45) - position.x > 0. && position.y - sourcePos.y * (1.-cos45) > 0. && sourcePos.y * (1.+cos45) - position.y > 0. && position.z - sourcePos.z * (1.-cos45) > 0. && sourcePos.z * (1.+cos45) - position.z > 0.) {
            isGood = field->getField(r / kpc, b);
        }
		else {
            //const Vector3d positionB = 2. * sourcePos + 2.*sourcePos.min() * sourceDis / sourceDis.getR() - position;
            //Vector3d sourceDisNew = sourceDis;
            Vector3d positionNew = position;
            //Vector3d positionB;
            while (positionNew.x - sourcePos.x * (1.-cos45) < 0. || sourcePos.x * (1.+cos45) - positionNew.x < 0.){
                if (positionNew.x - sourcePos.x * (1.-cos45) < 0.) positionNew.x = positionNew.x - 2*(positionNew.x - sourcePos.x * (1.-cos45));
                else positionNew.x = positionNew.x - 2*(positionNew.x - sourcePos.x * (1.+cos45));
            }
            while (positionNew.y - sourcePos.y * (1.-cos45) < 0. || sourcePos.y * (1.+cos45) - positionNew.y < 0.){
                if (positionNew.y - sourcePos.y * (1.-cos45) < 0.) positionNew.y = positionNew.y - 2*(positionNew.y - sourcePos.y * (1.-cos45));
                else positionNew.y = positionNew.y - 2*(positionNew.y - sourcePos.y * (1.+cos45));
            }
            while (positionNew.z - sourcePos.z * (1.-cos45) < 0. || sourcePos.z * (1.+cos45) - positionNew.z < 0.){
                if (positionNew.z - sourcePos.z * (1.-cos45) < 0.) positionNew.z = positionNew.z - 2*(positionNew.z - sourcePos.z * (1.-cos45));
                else positionNew.z = positionNew.z - 2*(positionNew.z - sourcePos.z * (1.+cos45));
            }
            /*
            while (sourceDisNew.x < 0. || sourceDisNew.y < 0. || sourceDisNew.z < 0. || sourceDisNew.x > sourcePos.min() || sourceDisNew.y > sourcePos.min() || sourceDisNew.z > sourcePos.min()){
                positionB = 2. * sourcePos + 2.*sourcePos.min() * sourceDisNew / sourceDisNew.getR() - positionNew;
                if (sourceDisNew.x < 0. || sourceDisNew.x > sourcePos.min()) positionNew.x = positionNew.x + 2*(sourcePos.min()-sourceDisNew.x); // positionNew.x = positionB.x; //
                if (sourceDisNew.y < 0. || sourceDisNew.y > sourcePos.min()) positionNew.y = positionNew.y + 2*(sourcePos.min()-sourceDisNew.y); //  positionNew.y = positionB.y;
                if (sourceDisNew.z < 0. || sourceDisNew.z > sourcePos.min()) positionNew.z = positionNew.z + 2*(sourcePos.min()-sourceDisNew.z); //  positionNew.z = positionB.z;
                //positionNew = positionB;
                sourceDisNew = positionNew - sourcePos;
            }
            * /
            /*
            const Vector3d NewsourceDis = positionB - sourcePos;
            if (NewsourceDis.getR()>sourcePos.min()) {
                isGood = field->getField(r / kpc, b);
            }
            */
            r = quimby::Vector3f(positionNew.x, positionNew.y, positionNew.z);
            isGood = field->getField(r / kpc, b);
		}
        
        if (!isGood) {
			std::ostringstream str;
			str << "ExtendedQuimbyMagneticField: invalid position at " << position;
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

} // namespace crpropa



#endif // CRPROPA_HAVE_QUIMBY
#endif // CRPROPA_QUIMBYMAGNETICFIELD_H
