/*******************************************************************************
 *  phong.cpp
 *******************************************************************************
 *  Copyright (c) 2013 Alexandre Kaspar <alexandre.kaspar@a3.epfl.ch>
 *  For Advanced Computer Graphics, at the LGG / EPFL
 *
 *        DO NOT REDISTRIBUTE
 ***********************************************/

#include <nori/acg.h>
#include <nori/bsdf.h>
#include <nori/frame.h>
#include "hemisampling.cpp"
#include <stdlib.h>
#include <time.h>
#include <cfloat>

NORI_NAMESPACE_BEGIN

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define GROUP_NUMBER 6
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

GROUP_NAMESPACE_BEGIN()

/**
 * \brief Phong BRDF model
 */
class Phong : public BSDF {
public:

    Phong(const PropertyList &propList) {
        m_Kd = propList.getColor("kd", Color3f(0.5f));
        m_Ks = propList.getColor("ks", Color3f(0.5f));
        m_exp = propList.getFloat("n", 20.0f);
        srand(time(NULL));
    }

    /// Reflection in local coordinates
    inline Vector3f reflect(const Vector3f &wi) const {
        return Vector3f(-wi.x(), -wi.y(), wi.z());
    }

    /// Evaluate the BRDF model

    Color3f eval(const BSDFQueryRecord &bRec) const {
        /* This is a smooth BRDF -- return zero if the measure
           is wrong, or when queried for illumination on the backside */
        if (bRec.measure != ESolidAngle
                || Frame::cosTheta(bRec.wi) <= 0
                || Frame::cosTheta(bRec.wo) <= 0) {
            return Color3f(0.0f);
        }

        // Based on http://www.cs.virginia.edu/~jdl/importance.doc‎
        float alpha = getAlpha(bRec);
        float kd = m_Kd.getLuminance();
        float ks = m_Ks.getLuminance();

        float result = kd/M_PI + ks * (m_exp+2.0f)/(2.0f*M_PI) * powf(cosf(alpha), m_exp);

        return Color3f(result);
    }

    /// Compute the density of \ref sample() wrt. solid angles

    float pdf(const BSDFQueryRecord &bRec) const {
        /* This is a smooth BRDF -- return zero if the measure
           is wrong, or when queried for illumination on the backside */
        if (bRec.measure != ESolidAngle
                || Frame::cosTheta(bRec.wi) <= 0
                || Frame::cosTheta(bRec.wo) <= 0)
            return 0.0f;

        float kd = m_Kd.getLuminance();
        float ks = m_Ks.getLuminance();
        float alpha = getAlpha(bRec);

        // compute both diffuse and specular pdfs
        float diffuseTerm = Frame::cosTheta(bRec.wo)/M_PI;
        float specularTerm = (m_exp + 1.0f) * powf(cosf(alpha), m_exp) / (2.0f * M_PI);

        // combine them using weighted average (kd + ks is probably one but... you never know)
        if (kd + ks == 0.0f) {
            return 0.0f;
        }
        float result = (kd * diffuseTerm + ks * specularTerm)/(kd + ks);
        // defensive programming
        if (result < 0.0f) {
            return 0.0f;
        }

        return result;
     }

    /// Draw a a sample from the BRDF model
    Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample_) const {
        Point2f sample(sample_);
        if (Frame::cosTheta(bRec.wi) <= 0)
            return Color3f(0.0f);

        float r = float(rand())/RAND_MAX;
        // 1. Select diffuse or specular
        float kd = m_Kd.getLuminance();
        float ks = m_Ks.getLuminance();
        float specSamplingWeight = ks / (ks + kd);
        bool useSpecular = true;
        if (r > specSamplingWeight) {
            useSpecular = false;
        }

        bRec.measure = ESolidAngle;
        bRec.eta = 1.0f; // no change in relative index

        float theta;
        float phi = 2.0f * M_PI * sample.y();
        if (useSpecular) {
            float exponent = 1.0f/(m_exp+1.0f);
            float alpha = acosf(powf(sample.x(), exponent));

            // This formula gives a vector in the reference frame of the specular direction
            theta = clamp(alpha, M_PI_2);
        } else {
            theta = acosf(sqrtf(sample.x()));
        }

        // Conversion to Cartesian coordinates
        float x = sinf(theta) * cosf(phi);
        float y = sinf(theta) * sinf(phi);
        float z = cosf(theta);
        bRec.wo = Vector3f(x, y, z);

        // Change of basis if the vector was defined in the reference frame of the specular direction
        if (useSpecular){
            convertToLocalFrame(bRec);
        }

        Color3f brdf = eval(bRec);
        float pdf_val = pdf(bRec);

        // If the vector is not probable at all, please ignore it, sir.
        if (pdf_val <= 0.0f) {
            return Color3f(0.0f);
        }
        return Color3f(brdf.getLuminance() / pdf_val) * Frame::cosTheta(bRec.wo);
    }

    /// Return a human-readable summary

    QString toString() const {
        return QString(
                "Phong[\n"
                "  Kd = %1\n"
                "  Ks = %2\n"
                "  n  = %3\n"
                "]").arg(m_Kd.toString()).arg(m_Ks.toString()).arg(m_exp);
    }

    Color3f getColor() const {
        return m_Kd;
    }

    EClassType getClassType() const {
        return EBSDF;
    }
private:
    Color3f m_Kd, m_Ks;
    float m_exp;

    inline float getAlpha(const BSDFQueryRecord &bRec) const {
        float cosTheta = Frame::cosTheta(bRec.wo);
        float cosThetaPlusAlpha = Frame::cosTheta(bRec.wi);

        float theta = acosf(cosTheta);
        float alpha = acosf(cosThetaPlusAlpha) - theta;

        return clamp(alpha, M_PI_2);
    }

    inline float clamp(const float f, const float limit) const {
        if (f > limit) {
            return limit;
        } else if (f < -limit) {
            return -limit;
        }
        return f;
    }

    inline void convertToLocalFrame(BSDFQueryRecord &bRec) const {
         Vector3f specularBase = bRec.wi;
         Frame specularRefFrame = Frame(specularBase);
         bRec.wo = specularRefFrame.toLocal(bRec.wo);
    }

};

GROUP_NAMESPACE_END

NORI_REGISTER_GROUP_CLASS(Phong, "phong");
NORI_NAMESPACE_END
