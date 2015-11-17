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

#include "Rigid_body_viewer.h"
#include <utils/gl.h>


//== IMPLEMENTATION ========================================================== 


Rigid_body_viewer::Rigid_body_viewer(const char* _title, int _width, int _height)
: Viewer_2D(_title, _width, _height)
{
    animate_             = false;
    time_step_           = 0.001;
    mass_                = 0.5;
    damping_linear_      = 0.1;
    damping_angular_     = 0.0001;
    spring_stiffness_    = 100.0;
    spring_damping_      = 5.0;

    mouse_spring_.active = false;
}


//-----------------------------------------------------------------------------


void Rigid_body_viewer::keyboard(int key, int x, int y)
{
    switch (key)
    {
    case '1':
    {
        std::vector<vec2> p;
        p.push_back( vec2(-0.6, -0.6) );
        p.push_back( vec2(-0.4, -0.6) );
        p.push_back( vec2(-0.4, -0.4) );
        p.push_back( vec2(-0.6, -0.4) );

        body_ = Rigid_body(p, mass_);
        body_.linear_velocity = vec2(5.0, 5.0);
        glutPostRedisplay();
        break;
    }
    case '2':
    {
        std::vector<vec2> p;
        p.push_back( vec2(-0.3, -0.1) );
        p.push_back( vec2(-0.1, -0.1) );
        p.push_back( vec2( 0.1, -0.1) );
        p.push_back( vec2( 0.3, -0.1) );
        p.push_back( vec2( 0.3,  0.1) );
        p.push_back( vec2( 0.1,  0.1) );
        p.push_back( vec2(-0.1,  0.1) );
        p.push_back( vec2(-0.3,  0.1) );

        body_ = Rigid_body(p, mass_);

        glutPostRedisplay();
        break;
    }
    case '3':
    {
        std::vector<vec2> p;
        p.push_back( vec2(-0.5,  0.1) );
        p.push_back( vec2(-0.5,  0.0) );
        p.push_back( vec2( 0.0,  0.0) );
        p.push_back( vec2( 0.0, -0.3) );
        p.push_back( vec2( 0.1, -0.3) );
        p.push_back( vec2( 0.1,  0.0) );
        p.push_back( vec2( 0.3,  0.0) );
        p.push_back( vec2( 0.3,  0.1) );

        body_ = Rigid_body(p, mass_);

        glutPostRedisplay();
        break;
    }
    // let parent class do the work
    default:
    {
        Viewer_2D::keyboard(key, x, y);
        break;
    }
    }
}


//-----------------------------------------------------------------------------


void Rigid_body_viewer:: mouse(int _button, int _state, int _x, int _y)
{
    // need points
    if (body_.points.empty())
        return;

    // mouse button release destroys current mouse spring
    if (_state == GLUT_UP)
    {
        mouse_spring_.active = false;
        return;
    }

    // mouse button press generates new mouse spring
    else if (_state == GLUT_DOWN)
    {
        // get point under mouse cursor
        vec2 p = pick(_x, _y);

        // find closest body point
        unsigned int i, imin;
        float dmin = FLT_MAX;
        for (i=0; i<body_.points.size(); ++i)
        {
            float d = distance(p, body_.points[i]);
            if (d < dmin)
            {
                dmin = d;
                imin = i;
            }
        }

        // setup the mouse spring
        mouse_spring_.active = true;
        mouse_spring_.particle_index = imin;
        mouse_spring_.mouse_position = p;
    }

    glutPostRedisplay();
}


//-----------------------------------------------------------------------------


void Rigid_body_viewer:: motion(int _x, int _y)
{
    if (mouse_spring_.active)
    {
        // update mouse position
        mouse_spring_.mouse_position = pick(_x, _y);
        glutPostRedisplay();
    }
}


//-----------------------------------------------------------------------------


void Rigid_body_viewer:: draw()
{
    // parent's status text
    Viewer_2D::draw();

    // draw walls
    glDisable(GL_LIGHTING);
    glLineWidth(1.0);
    glColor3f(0.5,0.5,0.5);
    glBegin(GL_LINE_STRIP);
    glVertex2f( -1.0,  1.0 );
    glVertex2f( -1.0, -1.0 );
    glVertex2f(  1.0, -1.0 );
    glVertex2f(  1.0,  1.0 );
    glVertex2f( -1.0,  1.0 );
    glEnd();

    // draw rigid body
    body_.draw();

    // draw mouse spring
    if (mouse_spring_.active)
    {
        glLineWidth(5.0);
        glColor3f(1,0,0);
        glBegin(GL_LINES);
        glVertex2fv( body_.points[ mouse_spring_.particle_index ].data() );
        glVertex2fv( mouse_spring_.mouse_position.data() );
        glEnd();
    }
}


//-----------------------------------------------------------------------------


void Rigid_body_viewer::compute_forces()
{ 
    /** \todo Compute all forces acting on the rigid body
     \li clear all forces
     \li add gravity
     \li add damping to linear and angular movement
     \li add the mouse spring force
     */
}


//-----------------------------------------------------------------------------


void Rigid_body_viewer::impulse_based_collisions()
{
    /** \todo Handle collisions based on impulses
     */
}


//-----------------------------------------------------------------------------


void Rigid_body_viewer::time_integration(float dt)
{
    // compute all forces
    compute_forces();

    /** \todo Implement explicit Euler time integration
     \li update position and orientation
     \li update linear and angular velocities
     \li call update_points() at the end to compute the new particle positions
     */
    // handle collisions
    impulse_based_collisions();
}


//=============================================================================
