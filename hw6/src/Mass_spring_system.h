//=============================================================================
//
//   Exercise code for the lecture
//   "Advanced Computer Graphics"
//
//   Adapted from Prof. Dr. Mario Botsch, Bielefeld University
//
//   Copyright (C) 2013 LGG, epfl
//
//   DO NOT REDISTRIBUTE
//=============================================================================

#ifndef MASS_SPRING_SYSTEM_H
#define MASS_SPRING_SYSTEM_H


//== INCLUDES =================================================================

#include "utils/vec2.h"
#include "utils/gl.h"
#include <vector>

//== CLASS DEFINITION =========================================================


/** \class Particle Mass_spring_system.h <01-mass_springs/Mass_spring_system.h>
 Class for a particle that stores position, velocity, force and other data.
 */
struct Particle
{
    /// constructor
    Particle(vec2 p, vec2 v, float m, bool l, int id_num)
    : position(p), velocity(v), force(0,0), mass(m), locked(l), acceleration(0,0), id (id_num)
    {}

    vec2  position; ///< position of the particle
    vec2  velocity; ///< velocity of the particle
    vec2  force;    ///< accumulated force acting on the particle
    float mass;     ///< mass of the particle
    bool  locked;   ///< is the particle locked?
    int id; ///< the id of the particle in the mass spring system

    vec2  position_t;   ///< used for Midpoint integration
    vec2  velocity_t;   ///< used for Midpoint integration
    vec2  acceleration; ///< used for Verlet integration
};



//== CLASS DEFINITION =========================================================


/** \class Spring Mass_spring_system.h <01-mass_springs/Mass_spring_system.h>
 Class for representing springs
 */
class Spring
{
public:
    /// construct with two particles. automatically computes rest length.
    Spring(Particle& p0, Particle& p1)
    : particle0(&p0), particle1(&p1)
    {
        rest_length = length();
    }

    /// get current length of the spring
    float length() const
    {
        return norm(particle0->position - particle1->position);
    }

    Particle*  particle0;   ///< pointer to first particle
    Particle*  particle1;   ///< pointer to second particle
    float      rest_length; ///< rest length
};



//== CLASS DEFINITION =========================================================


/** \class Triangle Mass_spring_system.h <01-mass_springs/Mass_spring_system.h>
 Class for storing triangles (for area preserving forces)
 */
struct Triangle
{
    /// construct with three particles. automatically computes rest area.
    Triangle(Particle& p0, Particle& p1, Particle& p2)
    : particle0(&p0), particle1(&p1), particle2(&p2)
    {
        rest_area = area();
    }

    /// compute current area
    float area() const
    {
        const vec2& p0 = particle0->position;
        const vec2& p1 = particle1->position;
        const vec2& p2 = particle2->position;
        return 0.5 * ((p1[0]-p0[0]) * (p2[1]-p0[1]) - (p2[0]-p0[0]) * (p1[1]-p0[1]));
    }

    Particle*  particle0; ///< pointer to first particle
    Particle*  particle1; ///< pointer to second particle
    Particle*  particle2; ///< pointer to third particle
    float      rest_area; ///< area in rest state
};



//== CLASS DEFINITION =========================================================


/** \class Mass_spring_system Mass_spring_system.h <01-mass_springs/Mass_spring_system.h>
 Class for managing a mass-spring system
 */
class Mass_spring_system
{
public:
    /// constructor
    Mass_spring_system() {}

    /// clear all particles, springs, and triangles
    void clear();

    /// add a particle
    void add_particle( vec2 position, vec2 velocity, float mass, bool locked );

    /// add a spring
    void add_spring(unsigned int i0, unsigned int i1);

    /// add a triangle
    void add_triangle(unsigned int i0, unsigned int i1, unsigned int i2);

    /// render the mass spring system
    void draw(float particle_radius, bool show_forces) const;

public:
    std::vector<Particle>  particles; ///< vector of all particles
    std::vector<Spring>    springs;   ///< vector of all springs
    std::vector<Triangle>  triangles; ///< vector of all triangles
};


//=============================================================================
#endif
//=============================================================================

