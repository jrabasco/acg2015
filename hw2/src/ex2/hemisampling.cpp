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

        // TODO implement cosine hemisphere sampling using a given uniform [0;1]^2 sample
        return Vector3f(0.0f, 0.0f, 1.0f); // this is wrong, replace it with your solution!
}

NORI_NAMESPACE_END
