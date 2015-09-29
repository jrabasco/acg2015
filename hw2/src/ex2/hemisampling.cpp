/*******************************************************************************
 *  hemisampling.cpp
 *******************************************************************************
 *  Copyright (c) 2013 Alexandre Kaspar <alexandre.kaspar@a3.epfl.ch>
 *  For Advanced Computer Graphics, at the LGG / EPFL
 * 
 *        DO NOT REDISTRIBUTE
 ***********************************************/

#include <nori/vector.h>

NORI_NAMESPACE_BEGIN

/**
 * \brief Cosine hemisphere sampling
 * 
 * \param sample a 2d uniform sample
 */
inline Vector3f squareToCosineHemisphere(const Point2f &sample) {
    // At first, we generate a random point on a sphere
    float randomRadius = sqrtf(sample.x());
    float randomAngle = 2 * M_PI * sample.y();

    // We apply pythagore knowing that sqrt(sample.x())^2 = sample.x()...
    float zCoordinate = sqrtf(1 - sample.x());

    return Vector3f(randomRadius * cosf(randomAngle),
                    randomRadius * sinf(randomAngle),
                    zCoordinate);
}

NORI_NAMESPACE_END
