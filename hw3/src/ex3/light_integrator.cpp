/*******************************************************************************
 *  light_integrator.cpp
 *******************************************************************************
 *  Copyright (c) 2013 Alexandre Kaspar <alexandre.kaspar@a3.epfl.ch>
 *  For Advanced Computer Graphics, at the LGG / EPFL
 * 
 *        DO NOT REDISTRIBUTE
 ***********************************************/

#include <nori/acg.h>
#include <nori/bsdf.h>
#include <nori/common.h>
#include <nori/integrator.h>
#include <nori/luminaire.h>
#include <nori/sampler.h>
#include <nori/scene.h>
#include <pthread.h>
#include <vector>

NORI_NAMESPACE_BEGIN

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define GROUP_NUMBER 6
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

GROUP_NAMESPACE_BEGIN()

/**
 * \brief Simple local illumination integrator
 * using light area sampling
 */
class LightIntegrator : public Integrator {
public:

        LightIntegrator(const PropertyList &propList) {
                Q_UNUSED(propList);
        }

        /// Return the mesh corresponding to a given luminaire
        inline const Mesh *getMesh(const Luminaire *lum) const {
                const Mesh *mesh = dynamic_cast<const Mesh *> (lum->getParent());
                if (!mesh) throw NoriException("Unhandled type of luminaire!");
                return mesh;
        }

        /**
         * \brief Directly sample the lights, providing a sample weighted by 1/pdf
         * where pdf is the probability of sampling that given sample
         * 
         * \param scene
         * the scene to work with
         * 
         * \param lRec
         * the luminaire information storage
         * 
         * \param _sample
         * the 2d uniform sample
         * 
         * \return the sampled light radiance including its geometric, visibility and pdf weights
         */
        inline Color3f sampleLights(const Scene *scene, LuminaireQueryRecord &lRec, const Point2f &_sample) const {
                Point2f sample(_sample);
                const std::vector<Luminaire *> &luminaires = scene->getLuminaires();

                if (luminaires.size() == 0)
                        throw NoriException("LightIntegrator::sampleLights(): No luminaires were defined!");

                // 1. Choose one luminaire at random
                int luminaireIdx = rand() % luminaires.size();
                lRec.luminaire = luminaires[luminaireIdx];

                // 2. Sample the position on the luminaire mesh
                // using Mesh::samplePosition(const Point2d &sample, Point3f &p, Normal3f &n)
                getMesh(lRec.luminaire)->samplePosition(sample, lRec.p, lRec.n);

                // 3. Compute geometry term G and visibility term on the luminaire's side (no cos(w) of the mesh side)
                // as well as the pdf of that point being found
                // use Mesh::pdf to get the probability of choosing the point in Mesh::samplePosition
                lRec = LuminaireQueryRecord(lRec.luminaire, lRec.ref, lRec.p, lRec.n);
                const Ray3f shadowRay(lRec.ref, lRec.d);

                Intersection its;
                const bool doesIntersect = scene->rayIntersect(shadowRay, its);
                if (!doesIntersect) {
                        if (scene->hasEnvLuminaire()) {
                               lRec = LuminaireQueryRecord(scene->getEnvLuminaire(), shadowRay);
                               return lRec.luminaire->eval(lRec);
                        } else {
                            return Color3f(0.0f);
                        }
                } else if (!its.mesh->isLuminaire()) {
                    return Color3f(0.0f);
                }

                const float cosThetaSecond = -lRec.d.dot(lRec.n);
                const bool doesLightPointsTowardsRef = cosThetaSecond > 0;
                if (!doesLightPointsTowardsRef) {
                    return 0.0f;
                }

                lRec.pdf = getMesh(lRec.luminaire)->pdf();

                // 4. Return radiance emitted from luminaire multiplied by the appropriate terms G, V ...
                if (lRec.pdf <= 0.0f || lRec.dist == 0.0f)
                    return 0.0f;


                const Color3f L_e = lRec.luminaire->eval(lRec);
                return ((L_e * cosThetaSecond) / (lRec.dist * lRec.dist)) / lRec.pdf;
        }

        /**
         * \brief Simple local illumination integration:
         * We cast a ray from the camera, intersects it with the first element
         * in the scene it goes through and from there we directly sample the
         * light's in the scene to compute direct lighting.
         */
        Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray_) const {
                Ray3f ray(ray_);

                /* Find the surface that is visible in the requested direction */
                Intersection its;
                if (!scene->rayIntersect(ray, its))
                        return Color3f(0.0f);

                const Mesh *mesh = its.mesh;
                const BSDF *bsdf = mesh->getBSDF();

                /* In case we hit a luminaire directly. */
                if (mesh->isLuminaire()) {
                    const Luminaire *luminaire = its.mesh->getLuminaire();
                    return luminaire->eval(LuminaireQueryRecord(luminaire, ray.o, its.p, its.shFrame.n));
                }


                LuminaireQueryRecord lRec;
                lRec.ref = its.p;
                const Point2f sample(sampler->next2D());
                Color3f sampledLight = sampleLights(scene, lRec, sample);

                float cosThetaPrime = Frame::cosTheta(ray.d);

                BSDFQueryRecord bRec(its.toLocal(-ray.d), its.toLocal(lRec.d), ESolidAngle);
                const Color3f f_r = bsdf->eval(bRec);

                return f_r * sampledLight * fabs(cosThetaPrime);
        }

        QString toString() const {
                return QString("LightIntegrator[]");
        }
};

GROUP_NAMESPACE_END

NORI_REGISTER_GROUP_CLASS(LightIntegrator, "light");
NORI_NAMESPACE_END
