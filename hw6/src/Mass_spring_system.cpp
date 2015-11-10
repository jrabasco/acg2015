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


//== INCLUDES =================================================================

#include "Mass_spring_system.h"


//== IMPLEMENTATION ========================================================== 


void Mass_spring_system::clear()
{
    particles.clear();
    springs.clear();
    triangles.clear();
}


//-----------------------------------------------------------------------------


void Mass_spring_system::add_particle(vec2 position, vec2 velocity, float mass, bool locked)
{
    particles.push_back( Particle(position, velocity, mass, locked, particles.size ()) );
}


//-----------------------------------------------------------------------------


void Mass_spring_system::add_spring(unsigned int i0, unsigned int i1)
{
    assert(i0 < particles.size());
    assert(i1 < particles.size());
    springs.push_back( Spring(particles[i0], particles[i1]) );
}


//-----------------------------------------------------------------------------


void Mass_spring_system::add_triangle(unsigned int i0, unsigned int i1, unsigned int i2)
{
    assert(i0 < particles.size());
    assert(i1 < particles.size());
    assert(i2 < particles.size());
    triangles.push_back( Triangle(particles[i0], particles[i1], particles[i2]) );
}


//-----------------------------------------------------------------------------


void Mass_spring_system::draw(float particle_radius, bool show_forces) const
{
    // draw particles
    glEnable(GL_LIGHTING);
    glShadeModel(GL_SMOOTH);
    glEnable(GL_COLOR_MATERIAL);
    for (unsigned int i=0; i<particles.size(); ++i)
    {
        const Particle& p = particles[i];
        if (p.locked)
            glColor3f(0.6, 0.0, 0.0);
        else
            glColor3f(0.0, 0.6, 0.0);

        glPushMatrix();
        glTranslatef( p.position[0], p.position[1], 0.0 );
        glutSolidSphere( particle_radius, 20, 20 );
        glPopMatrix();

        if (show_forces)
        {
            glLineWidth(1.0);
            glBegin(GL_LINES);
            glVertex2f(p.position[0], p.position[1]);
            glVertex2f(p.position[0] + 0.05 * p.force[0], p.position[1] + 0.05 * p.force[1]);
            glEnd();
        }
    }
    glDisable(GL_COLOR_MATERIAL);

    // draw springs
    glDisable(GL_LIGHTING);
    glLineWidth(5.0);
    glColor3f(0,0,0);
    glBegin(GL_LINES);
    for (unsigned int i=0; i<springs.size(); ++i)
    {
        glVertex2fv( springs[i].particle0->position.data() );
        glVertex2fv( springs[i].particle1->position.data() );
    }
    glEnd();

    // draw area constraints
    glDisable(GL_LIGHTING);
    glColor3f(0.8, 1.0, 0.8);
    glBegin(GL_TRIANGLES);
    for (unsigned int i=0; i<triangles.size(); ++i)
    {
        glVertex2fv( triangles[i].particle0->position.data() );
        glVertex2fv( triangles[i].particle1->position.data() );
        glVertex2fv( triangles[i].particle2->position.data() );
    }
    glEnd();
}

//=============================================================================
