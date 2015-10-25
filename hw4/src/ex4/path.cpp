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
                float Q = 0.2f; // for Russian Roulette, set to 1 to disable

                // Variables used for recursion (though it's not technically recursion)
                Ray3f ray(_ray);
                Color3f result(0.0f), throughput(1.0f);

                for(int n = 1; ; ++n) {
                    Intersection its;
                    // STEP 1
                    // Intersect the ray with the scene. Return environment luminaire if no hit.
                    if(!scene->rayIntersect(ray, its)) {
                        if(scene->hasEnvLuminaire()) {
                            LuminaireQueryRecord envRec(scene->getEnvLuminaire(), ray);
                            result += throughput * envRec.luminaire->eval(envRec);
                        }

                        break;
                    }

                    // STEP 2
                    // In case w_i intersects with a luminaire, we take it into account only in case this is the first bounce
                    // Otherwise, we have already taken it into account during the last direct lighting phase.
                    if(its.mesh->isLuminaire() && n == 1) {
                        LuminaireQueryRecord lRec(its.mesh->getLuminaire(), ray.o, its.p, its.shFrame.n);
                        result += throughput * lRec.luminaire->eval(lRec);
                    }

                    // STEP 3
                    // Compute direct lighting on the intersection point by sampling the lights
                    LuminaireQueryRecord lRec(its.p);
                    Vector3f w_i = its.toLocal(-ray.d);
                    Color3f directColor = sampleLights(scene, lRec, sampler->next2D());
                    BSDFQueryRecord bsdfRec(w_i, its.toLocal(lRec.d), ESolidAngle);
                    result += throughput * its.mesh->getBSDF()->eval(bsdfRec) * directColor * std::abs(Frame::cosTheta(bsdfRec.wo));


                    // STEP 4
                    // Sample the BSDF to get compute the next direction of the path
                    BSDFQueryRecord bsdfSampleRec(w_i);
                    Color3f sampledColor = its.mesh->getBSDF()->sample(bsdfSampleRec, sampler->next2D());

                    // In case we sampled 0x000000 or to avoid division by 0
                    if ((sampledColor.array() == 0).all() || its.mesh->pdf() == 0) {
                        break;
                    }

                    throughput *= sampledColor;
                    ray = Ray3f(its.p, its.shFrame.toWorld(bsdfSampleRec.wo));

                    // Step 5. Apply Russian Roulette after 2 main bounces.
                    if(n >= SAMPLE_DEPTH) {
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

private:
    static const int SAMPLE_DEPTH = 2;
};

GROUP_NAMESPACE_END

NORI_REGISTER_GROUP_CLASS(PathTracer, "path");
NORI_NAMESPACE_END
