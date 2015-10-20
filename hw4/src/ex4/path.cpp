/*******************************************************************************
 *  path.cpp
 *******************************************************************************
 *  Copyright (c) 2013 Alexandre Kaspar <alexandre.kaspar@a3.epfl.ch>
 *  For Advanced Computer Graphics, at the LGG / EPFL
 *
 *        DO NOT REDISTRIBUTE
 ***********************************************/

#include <nori/acg.h>
#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/luminaire.h>
#include <nori/bsdf.h>
#include <nori/medium.h>
#include <nori/phase.h>
#include <nori/scene.h>

NORI_NAMESPACE_BEGIN

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define GROUP_NUMBER 6
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

GROUP_NAMESPACE_BEGIN()

/**
 * Simple path tracer implementation
 */
class PathTracer : public Integrator {
public:

        PathTracer(const PropertyList &) {
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
                int index = std::min((int) (luminaires.size() * sample.x()), (int) luminaires.size() - 1);
                sample.x() = luminaires.size() * sample.x() - index; // process sample to be Unif[0;1] again

                // 2. Sample the position on the luminaire mesh
                // using Mesh::samplePosition(const Point2d &sample, Point3f &p, Normal3f &n)
                lRec.luminaire = luminaires[index];
                const Mesh *mesh = getMesh(lRec.luminaire);
                mesh->samplePosition(sample, lRec.p, lRec.n);
                lRec.d = lRec.p - lRec.ref;

                // 3. Compute distance between the two points (from first mesh, to luminaire mesh)
                float dist2 = lRec.d.squaredNorm();
                lRec.dist = std::sqrt(dist2);
                lRec.d /= lRec.dist;

                // 4. Correct side of luminaire
                // /!\ if on the wrong side, then we get no contribution!
                float dp = -lRec.n.dot(lRec.d);
                lRec.pdf = dp > 0 ? mesh->pdf() * dist2 / dp : 0.0f;

                if (dp > 0) {
                        // 5. Check the visibility
                        if (scene->rayIntersect(Ray3f(lRec.ref, lRec.d, Epsilon, lRec.dist * (1 - 1e-4f))))
                                return Color3f(0.0f);
                        // 6. Geometry term on luminaire's side
                        // Visiblity + Geometric term on the luminaire's side 
                        //      G(x, x', w, w') = ( cos(w) cos(w') ) / ||x - x'||^2
                        float G_lum = dp / dist2;

                        // 7. Radiance from luminaire
                        Color3f value = lRec.luminaire->getColor();

                        return value * G_lum * luminaires.size() / mesh->pdf();
                } else {
                        // wrong side of luminaire!
                        return Color3f(0.0f);
                }
        }

        Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &_ray) const {
                float Q = 1; // for Russian Roulette, set to 1 to disable

                // Variables used for recursion (though it's not technically recursion)
                Ray3f ray(_ray);
                Color3f result(0.0f), throughput(1.0f);

                for(int n = 0; ; ++n) {
                    Intersection its;
                    // Step 1: Intersect the ray with the scene. Return environment luminaire if no hit.
                    if(!scene->rayIntersect(ray, its)) {
                        // TODO: Write (obvious since it's in the assignment) justification for that in the report
                        // TODO: Why is this false sometimes?
                        if(scene->hasEnvLuminaire()) {
                            LuminaireQueryRecord envRec(scene->getEnvLuminaire(), ray);
                            result += throughput * envRec.luminaire->eval(envRec);
                        }

                        break;

                    }

                    // Step 2: Check if the ray hit a light source.
                    if(its.mesh->isLuminaire()) {
                        // TODO: Write (relatively obvious IMHO) justification for that in the report
                        LuminaireQueryRecord meshRec(its.mesh->getLuminaire(), ray.o, its.p, its.shFrame.n);
                        // TODO: This line causes a metric tone of light noise in the resulting image, something's wrong
                        result += throughput * meshRec.luminaire->eval(meshRec);
                        break;
                    }

                    // Step 3: Direct illumination sampling. (i.e. n = 0)
                    LuminaireQueryRecord lRec(its.p);
                    Color3f directColor = sampleLights(scene, lRec, sampler->next2D());
                    Vector3f w_i = its.toLocal(-ray.d); // Invert the ray since we want w_i going *from* the point
                    // N.B.: BSDFQueryRecord works in local coords while this function works in global coords,
                    //       so we must convert every time.
                    BSDFQueryRecord bsdfRec(w_i, its.toLocal(lRec.d), ESolidAngle);
                    // Create a ray from the light-intersecting info
                    Ray3f transmittanceRay(lRec.ref, lRec.d, 0, lRec.dist);
                    // TODO: I think this is missing the magical G thingy
                    result += throughput * scene->evalTransmittance(transmittanceRay, sampler) * directColor * its.mesh->getBSDF()->eval(bsdfRec);

                    // Step 4: Recursively sample indirect illumination (i.e. n > 0)
                    BSDFQueryRecord bsdfSampleRec(w_i);
                    Color3f sampledColor = its.mesh->getBSDF()->sample(bsdfSampleRec, sampler->next2D());
                    if(sampledColor.getLuminance() == 0) {
                        // bail out early if sampledColor is black, no point in continuing with throughput == 0
                        // without this check it gets really really slow
                        break;
                    }
                    throughput *= sampledColor;
                    ray = Ray3f(its.p, its.shFrame.toWorld(bsdfSampleRec.wo));

                    // Step 5. Apply Russian Roulette after 2 main bounces.
                    if(n > 2) {
                        float random = sampler->next1D();
                        if(random > Q) {
                            throughput /= (1 - Q);
                        } else {
                            break;
                        }
                    }
                }
                return result;
        }

        QString toString() const {
                return "PathTracer[]";
        }
};

GROUP_NAMESPACE_END

NORI_REGISTER_GROUP_CLASS(PathTracer, "path");
NORI_NAMESPACE_END
