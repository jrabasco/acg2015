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

#include "Mass_spring_viewer.h"
#include "utils/gl_text.h"
#include <sstream>
#include <fstream>
#include <string>

bool equilibrium = false;
bool output = false;

//== IMPLEMENTATION ==========================================================

float planes[4][3] = {
    {  0.0,  1.0, 1.0 },
    {  0.0, -1.0, 1.0 },
    {  1.0,  0.0, 1.0 },
    { -1.0,  0.0, 1.0 }
};


Mass_spring_viewer::
Mass_spring_viewer(const char* _title, int _width, int _height)
: Viewer_2D(_title, _width, _height)
{
    integration_         = Euler;
    collisions_          = Force_based;
    external_force_      = None;
    animate_             = false;
    area_forces_         = false;
    show_forces_         = false;

    time_step_           = 0.001;
    particle_radius_     = 0.03;

    particle_mass_       = 0.1;
    damping_             = 0.1;
    collision_stiffness_ = 1000.0;
    spring_stiffness_    = 1000.0;
    spring_damping_      = 1.0;
    area_stiffness_      = 100000.0;

    mouse_spring_.active = false;
}


//-----------------------------------------------------------------------------


void Mass_spring_viewer::keyboard(int key, int x, int y)
{
    switch (key)
    {
        // setup problem 1
        case '1':
        {
            body_.clear();
            body_.add_particle( vec2(-0.5, -0.5), vec2(14.0, -2.0), particle_mass_, false );
            glutPostRedisplay();
            break;
        }

        // setup problem 2
        case '2':
        {
            body_.clear();
            for (int i=0; i<100; ++i)
            {
                body_.add_particle( vec2(0.9* cos(i/50.0*M_PI), 0.9*sin(i/50.0*M_PI)), vec2(-sin(i/50.0*M_PI), cos(i/50.0*M_PI)), particle_mass_, false );
            }

            glutPostRedisplay();
            break;
        }

        // setup problem 3
        case '3':
        {
            body_.clear();

            for (int i=0; i<10; ++i)
            {
                body_.add_particle( vec2(i*0.1, 0.8), vec2(0.0, 0.0), particle_mass_, i==0 );
            }

            for (unsigned int i=0; i<9; ++i)
            {
                body_.add_spring(i, i+1);
            }

            glutPostRedisplay();
            break;
        }

        // setup problem 4
        case '4':
        {
            body_.clear();

            body_.add_particle( vec2(-0.1, 0.7), vec2(0.0, 0.0), particle_mass_, false );
            body_.add_particle( vec2( 0.0, 0.6), vec2(0.0, 0.0), particle_mass_, false );
            body_.add_particle( vec2( 0.1, 0.7), vec2(0.0, 0.0), particle_mass_, false );

            body_.add_spring(0, 1);
            body_.add_spring(0, 2);
            body_.add_spring(1, 2);

            body_.add_triangle(0, 1, 2);

            glutPostRedisplay();
            break;
        }

        // setup problem 5
        case '5':
        {
            body_.clear();

            for (int i=0; i<8; ++i)
            {
                body_.add_particle( vec2(-0.5+0.2*cos(0.25*i*M_PI), -0.5+0.2*sin(0.25*i*M_PI)), vec2(5.0, 5.0), particle_mass_, false );
            }

            body_.add_particle( vec2(-0.5, -0.5), vec2(5.0, 5.0), particle_mass_, false );

            for (unsigned int i=0; i<8; ++i)
            {
                body_.add_spring( i, (i+1) % 8 );
                body_.add_spring(i, 8);
                body_.add_triangle(i, (i+1)%8, 8);
            }

            glutPostRedisplay();
            break;
        }

        // switch between time integration methods
        case 'i':
        {
            switch (integration_)
            {
                case Euler:    integration_ = Midpoint; break;
                case Midpoint: integration_ = Verlet;   break;
                case Verlet:   integration_ = Implicit;    break;
                case Implicit: integration_ = Euler; break;
            }
            glutPostRedisplay();
            break;
        }

        // switch between center and gravitation force
        case 'f':
        {
            switch (external_force_)
            {
                case None:        external_force_ = Center; break;
                case Center:      external_force_ = Gravitation; break;
                case Gravitation: external_force_ = None;      break;
            }
            glutPostRedisplay();
            break;
        }

        // switch between force-based and impulse-based collisions
        case 'c':
        {
            switch (collisions_)
            {
                case Force_based:   collisions_ = Impulse_based; break;
                case Impulse_based: collisions_ = Force_based;   break;
            }
            glutPostRedisplay();
            break;
        }


        // toggle area forces on/off
        case 'a':
        {
            area_forces_ = !area_forces_;
            glutPostRedisplay();
            break;
        }


        // visualization of particle forces on/off
        case 'v':
        {
            show_forces_ = !show_forces_;
            glutPostRedisplay();
            break;
        }

        case 'e':
        {
            equilibrium = !equilibrium;
            break;
        }

        case 'o':
        {
            output = !output;
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


void Mass_spring_viewer::draw()
{
    // parent's status text
    Viewer_2D::draw();


    // draw some status text
    glDisable(GL_LIGHTING);
    glColor3f(1,0,0);
    std::ostringstream oss;

    oss.str("");
    oss << "Integration: ";
    switch (integration_)
    {
        case Euler:    oss << "Euler";    break;
        case Midpoint: oss << "Midpoint"; break;
        case Verlet:   oss << "Verlet";   break;
        case Implicit: oss << "Implicit"; break;
    }
    glText(20, height_-40, oss.str());

    oss.str("");
    oss << "#Particles: " << body_.particles.size();
    glText(20, height_-60, oss.str());

    oss.str("");
    oss << "#Springs: " << body_.springs.size();
    glText(20, height_-80, oss.str());

    oss.str("");
    oss << "#Triangles: " << body_.triangles.size();
    glText(20, height_-100, oss.str());

    oss.str("");
    oss << "Area Forces: " << (area_forces_ ? "on" : "off");
    glText(20, height_-120, oss.str());

    oss.str("");
    oss << "Collisions: " << (collisions_ == Force_based ? "force" : "impulse");
    glText(20, height_-140, oss.str());

    oss.str("");
    oss << "External force: ";
    switch (external_force_)
    {
        case None:        oss << "None";        break;
        case Center:      oss << "Center";      break;
        case Gravitation: oss << "Gravitation"; break;
    }
    glText(20, height_-160, oss.str());

    oss.str("");
    oss << "Visualize forces: " << (show_forces_ ? "on" : "off");
    glText(20, height_-180, oss.str());


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


    // draw mouse spring
    if (mouse_spring_.active)
    {
        glDisable(GL_LIGHTING);
        glLineWidth(5.0);
        glColor3f(1,0,0);
        glBegin(GL_LINES);
        glVertex2fv( body_.particles[mouse_spring_.particle_index].position.data() );
        glVertex2fv( mouse_spring_.mouse_position.data() );
        glEnd();
    }


    // draw particles, springs, triangles
    body_.draw(particle_radius_, show_forces_);
}


//-----------------------------------------------------------------------------


void Mass_spring_viewer::mouse(int _button, int _state, int _x, int _y)
{
    // need particles to do interaction
    if (!body_.particles.empty())
    {
        // mouse button release destroys current mouse spring
        if (_state == GLUT_UP)
        {
            mouse_spring_.active = false;
        }

        // mouse button press generates new mouse spring
        else if (_state == GLUT_DOWN)
        {
            // get point under mouse cursor
            vec2 p = pick(_x, _y);

            // find closest particle
            int   pidx = -1;
            float dmin = FLT_MAX;
            for (unsigned int i=0; i<body_.particles.size(); ++i)
            {
                float d = norm(p - body_.particles[i].position);
                if (d < dmin)
                {
                    dmin = d;
                    pidx = i;
                }
            }

            // construct mouse spring
            mouse_spring_.mouse_position = p;
            mouse_spring_.particle_index = pidx;
            mouse_spring_.active = true;
        }
    }

    glutPostRedisplay();
}


//-----------------------------------------------------------------------------


void Mass_spring_viewer::motion(int _x, int _y)
{
    if (mouse_spring_.active)
    {
        // update mouse positions
        mouse_spring_.mouse_position = pick(_x, _y);
    }

    glutPostRedisplay();
}


//-----------------------------------------------------------------------------


void Mass_spring_viewer::time_integration(float dt)
{
    std::vector<Particle>& particles = body_.particles;
    switch (integration_)
    {
        case Euler:
        {
            /** \todo (Part 1) Implement Euler integration scheme
             \li The Particle class has variables position_t and velocity_t to store current values
             \li Hint: compute_forces() computes all forces for the current positions and velocities.
             */
            compute_forces();
            for (std::vector<Particle>::iterator particle = particles.begin(); particle != particles.end(); ++particle) {
                if (!particle->locked) {
                    particle->acceleration = particle->force / particle->mass;
                    particle->velocity += dt * particle->acceleration;
                    particle->position += dt * particle->velocity;
                }
            }

            break;
        }

        case Midpoint:
        {
            /** \todo (Part 2) Implement the Midpoint time integration scheme
             \li The Particle class has variables position_t and velocity_t to store current values
             \li Hint: compute_forces() computes all forces for the current positions and velocities.
             */
            compute_forces();
            for (std::vector<Particle>::iterator particle = particles.begin(); particle != particles.end(); ++particle) {
                if (!particle->locked) {
                    particle->position_t = particle->position;
                    particle->velocity_t = particle->velocity;
                    particle->acceleration = particle->force / particle->mass;
                    particle->position += (dt / 2.0) * particle->velocity;
                    particle->velocity += (dt / 2.0) * particle->acceleration;
                }
            }
            compute_forces();
            for (std::vector<Particle>::iterator particle = particles.begin(); particle != particles.end(); ++particle) {
                if (!particle->locked) {
                    particle->acceleration = particle->force / particle->mass;
                    particle->position = particle->position_t + dt * particle->velocity;
                    particle->velocity = particle->velocity_t + dt * particle->acceleration;
                }
            }

            break;
        }


        case Verlet:
        {
            /** \todo (Part 2) Implement the Verlet time integration scheme
             \li The Particle class has a variable acceleration to remember the previous values
             \li Hint: compute_forces() computes all forces for the current positions and velocities.
             */
            compute_forces();
            for (std::vector<Particle>::iterator particle = particles.begin(); particle != particles.end(); ++particle) {
                if (!particle->locked) {
                    particle->acceleration = particle->force / particle->mass;
                    vec2 x = particle->position;
                    particle->position += dt * particle->velocity + dt * dt / 2.0 * particle->acceleration;
                    particle->velocity_t = particle->velocity;
                    particle->velocity = (particle->position - x) / dt;
                }
            }
            compute_forces();
            for (std::vector<Particle>::iterator particle = particles.begin(); particle != particles.end(); ++particle) {
                if (!particle->locked) {
                    vec2 a = particle->acceleration;
                    particle->acceleration = particle->force / particle->mass;
                    particle->velocity = particle->velocity_t + dt * (a + particle->acceleration) / 2.0;
                }
            }

            break;
        }

        case Implicit:
        {
            /// The usual force computation method is called, and then the jacobian matrix dF/dx is calculated
            compute_forces ();
            compute_jacobians ();

            /// Finally the linear system is composed and solved
            solver_.solve (dt, particle_mass_, body_.particles);

            break;
        }
    }


    // impulse-based collision handling
    if (collisions_ == Impulse_based)
    {
        impulse_based_collisions();
    }


    glutPostRedisplay();
}


//-----------------------------------------------------------------------------


void
Mass_spring_viewer::compute_forces()
{
    // clear forces
        for (unsigned int i=0; i<body_.particles.size(); ++i)
            body_.particles[i].force = vec2(0,0);


        std::vector<Particle>& particles = body_.particles;
        for (std::vector<Particle>::iterator particle = particles.begin(); particle != particles.end(); ++particle) {
            if (!particle->locked) {

                /** \todo (Part 1) Implement center force
                 */
                if (external_force_ == Center)
                {
                    static const float c = 20.0;
                    particle->force -= c * particle->position;
                }


                /** \todo (Part 1) Implement damping force
                 \li The damping coefficient is given as member damping_
                 */
                 if (integration_ != Implicit)
                    particle->force -= damping_*particle->velocity;


                /** \todo (Part 1) Implement gravitation force
                 \li Particle mass available as particle_mass_
                 */
                if (external_force_ == Gravitation)
                {
                    particle->force -= vec2(0, 9.81 * particle->mass);
                }


                /** \todo (Part 1) Implement force based boundary collisions
                 \li Collision coefficient given as collision_stiffness_
                 */
                // collision forces
                if (collisions_ == Force_based)
                {

                    for (int i = 0; i < 4; ++i) {
                        float A = planes[i][0], B = planes[i][1], C = -planes[i][2];
                        vec2 n(A, B);
                        vec2 p;
                        if (A != 0.0) {
                            p = vec2(C / A, 0.0);
                        } else if (B != 0.0) {
                            p = vec2(0.0, C / B);
                        }

                        float d = fabs(A * particle->position[0] + B * particle->position[1] - C)/sqrtf(A * A + B * B);

                        if (dot(p - particle->position, n) >= 0 || d <= particle_radius_) {
                            particle->force += collision_stiffness_ * (d + particle_radius_) * n;
                        }

                    }


                }
            }
        }


        /** \todo (Part 1) Compute force of the interactive mouse spring
         \li Required coefficients are given as spring_stiffness_ and spring_damping_
         */
        if (mouse_spring_.active)
        {
            Particle& p0 = body_.particles[ mouse_spring_.particle_index ];

            if (!p0.locked) {
                vec2 pos0 = p0.position;
                vec2 pos1 = mouse_spring_.mouse_position;

                vec2 direction = pos0 - pos1;
                float stretchLength = norm(direction);
                vec2 norm_direction = direction/stretchLength;

                p0.force -= (spring_stiffness_ * stretchLength  + spring_damping_ * dot(p0.velocity, norm_direction)) * norm_direction;
            }
        }


        /** \todo (Part 1) Compute spring forces
         \li Required information about springs in the scene are found in body_.springs
         \li Required coefficients are given as spring_stiffness_ and spring_damping_
         */

        std::vector<Spring>& springs = body_.springs;

        for (std::vector<Spring>::iterator spring = springs.begin(); spring != springs.end(); ++spring) {
            Particle* p0 = spring->particle0;
            Particle* p1 = spring->particle1;
            if (!p0->locked || !p1->locked) {
                float ks = spring_stiffness_;
                float kd = spring_damping_;
                float L = spring->rest_length;
                float xx0 = spring->length();
                vec2 norm_direction = (p0->position - p1->position) / xx0;
                vec2 v_diff = p0->velocity - p1->velocity;
                vec2 F0 = -(ks * (xx0 - L) + kd * dot(v_diff, norm_direction)) * norm_direction;
                if (!p0->locked) {
                    p0->force += F0;
                }
                if (!p1->locked) {
                    p1->force -= F0;
                }
            }
        }

    /** \todo (Part 2) Compute more forces in part 2 of the exercise: triangle-area forces, binding forces, etc.
     */
    if (area_forces_)
    {
        std::vector<Triangle>& triangles = body_.triangles;
        for (std::vector<Triangle>::iterator triangle = triangles.begin(); triangle != triangles.end(); ++triangle) {
            const vec2& p0 = triangle->particle0->position;
            const vec2& p1 = triangle->particle1->position;
            const vec2& p2 = triangle->particle2->position;
            const float coeff = area_stiffness_ * (triangle->area() - triangle->rest_area) / 2.0;
            vec2 dArea0(p1[1] - p2[1], p2[0] - p1[0]);
            vec2 dArea1(p2[1] - p0[1], p0[0] - p2[0]);
            vec2 dArea2(p0[1] - p1[1], p1[0] - p0[0]);
            triangle->particle0->force -= coeff * dArea0;
            triangle->particle1->force -= coeff * dArea1;
            triangle->particle2->force -= coeff * dArea2;
        }
    }

    if (equilibrium)
    {
        for (std::vector<Particle>::iterator particle = particles.begin(); particle != particles.end(); ++particle) {
            if (!particle->locked) {
                for (std::vector<Particle>::iterator p2 = particles.begin(); p2 != particles.end(); ++p2) {
                    if (p2 != particle) {
                        vec2 dir = particle->position - p2->position;
                        float dist = norm(dir);
                        dir = dir/dist;
                        particle->force += 0.01 * (dir/(dist * dist));

                    }
                }
                for (int i = 0; i < 4; ++i) {
                    float A = planes[i][0], B = planes[i][1], C = planes[i][2];
                    vec2 n(A, B);
                    float d = A*particle->position[0] + B*particle->position[1] + C;
                    particle->force += 0.1 * n / (d*d);
                }
            }
        }
    }

    if (output) {
    float kinetic_energy = 0.0f;
    float potential_energy = 0.0f;
    std::ofstream out;
    switch (integration_)
    {
        case Euler:
        {
            out.open("euler.csv", std::ios::app);
            break;
        }
        case Midpoint:
        {
            out.open("midpoint.csv", std::ios::app);
            break;
        }
        case Verlet:
        {
            out.open("verlet.csv", std::ios::app);
            break;
        }
        case Implicit:
        {
            out.open("implicit.csv", std::ios::app);
            break;
        }
    }

    for (std::vector<Particle>::iterator particle = particles.begin(); particle != particles.end(); ++particle) {
        float speed = norm(particle->velocity);
        kinetic_energy += particle->mass * speed * speed / 2.0f;
        potential_energy += particle->mass * (particle->position[1] + 1.0f) * 9.81;
    }

    out << kinetic_energy << "\n";
    out.close();
    }
}


//-----------------------------------------------------------------------------


void Mass_spring_viewer::impulse_based_collisions()
{
    /** \todo (Part 2) Handle collisions based on impulses
     */

    float epsilon = 0.9;
    std::vector<Particle>& particles = body_.particles;
    for (std::vector<Particle>::iterator particle = particles.begin(); particle != particles.end(); ++particle) {
        for (int i = 0; i < 4; ++i) {
            float A = planes[i][0], B = planes[i][1], C = planes[i][2];
            vec2 n(A, B);
            float d = A*particle->position[0] + B*particle->position[1] + C;
            float vn = dot(n, particle->velocity);
            if (d <= particle_radius_ && vn < 0.0) {
                particle->velocity -= (1 + epsilon) * n * vn;
            }
        }
    }
}
//=============================================================================

void Mass_spring_viewer::compute_jacobians ()
{
  /// Clear the solver matrices
  solver_.clear ();

  std::vector<Particle>& particles = body_.particles;

  // Center forces
  if (external_force_ == Center) {
    static const float c = 20.0;
    for (std::vector<Particle>::iterator particle = particles.begin(); particle != particles.end(); ++particle) {
        if (!particle->locked) {
            int i = particle->id;
            solver_.addElementToJacobian(2 * i, 2 * i, -c);
            solver_.addElementToJacobian(2 * i + 1, 2 * i + 1, -c);
        }
    }
  }

  // Collisions
  if (collisions_ == Force_based) {
    for (std::vector<Particle>::iterator particle = particles.begin(); particle != particles.end(); ++particle) {
        if (particle->locked)
            continue;

        for (int i = 0; i < 4; ++i) {
            float A = planes[i][0], B = planes[i][1], C = planes[i][2];
            vec2 n(A, B);
            float d = A*particle->position[0] + B*particle->position[1] + C;
            if (d <= particle_radius_) {
                float sgn = (d < 0) ? -1.0f : 1.0f;
                float coeff = sgn * collision_stiffness_;
                float jacob_idx = 2 * particle->id;
                solver_.addElementToJacobian(jacob_idx,     jacob_idx,     coeff * A * n[0]);
                solver_.addElementToJacobian(jacob_idx,     jacob_idx + 1, coeff * B * n[0]);
                solver_.addElementToJacobian(jacob_idx + 1, jacob_idx,     coeff * A * n[1]);
                solver_.addElementToJacobian(jacob_idx + 1, jacob_idx + 1, coeff * B * n[1]);
            }
        }

    }
  }

  // Springs
  std::vector<Spring>& springs = body_.springs;

  for (std::vector<Spring>::iterator spring = springs.begin(); spring != springs.end(); ++spring) {
    Particle * const pi = spring->particle0;
    Particle * const pj = spring->particle1;

    const int i = pi->id;
    const int j = pj->id;
    assert(i < particles.size());
    assert(j < particles.size());

    const vec2 xi = pi->position;
    const vec2 xj = pj->position;
    const vec2 diff = xi - xj;
    const float dist = spring->length();
    const float distSq = dist * dist;

    const vec2 I2[2] = {
        vec2(1.0f, 0.0f),
        vec2(0.0f, 1.0f)
    };

    const vec2 diffDiffT[2] = {
        vec2(diff[0] * diff[0], diff[0] * diff[1]),
        vec2(diff[0] * diff[1], diff[1] * diff[1])
    };

    vec2 dFi_dxi[2] = {
        -spring_stiffness_ * (I2[0] - spring->rest_length * (I2[0] - diffDiffT[0] / distSq) / dist) - spring_damping_ * diffDiffT[0] / (distSq * time_step_),
        -spring_stiffness_ * (I2[1] - spring->rest_length * (I2[1] - diffDiffT[1] / distSq) / dist) - spring_damping_ * diffDiffT[1] / (distSq * time_step_)
    };

    vec2 dFi_dxj[2] = {
        -dFi_dxi[0],
        -dFi_dxi[1]
    };

    if (!pi->locked) {
      solver_.addElementToJacobian(2 * i,     2 * i,      dFi_dxi[0][0]);
      solver_.addElementToJacobian(2 * i,     2 * i + 1,  dFi_dxi[0][1]);
      solver_.addElementToJacobian(2 * i + 1, 2 * i,      dFi_dxi[1][0]);
      solver_.addElementToJacobian(2 * i + 1, 2 * i + 1,  dFi_dxi[1][1]);

      solver_.addElementToJacobian(2 * i,     2 * j,      dFi_dxj[0][0]);
      solver_.addElementToJacobian(2 * i,     2 * j + 1,  dFi_dxj[0][1]);
      solver_.addElementToJacobian(2 * i + 1, 2 * j,      dFi_dxj[1][0]);
      solver_.addElementToJacobian(2 * i + 1, 2 * j + 1,  dFi_dxj[1][1]);
    }

    if (!pj->locked) {
      solver_.addElementToJacobian(2 * j,     2 * i,     -dFi_dxi[0][0]);
      solver_.addElementToJacobian(2 * j,     2 * i + 1, -dFi_dxi[0][1]);
      solver_.addElementToJacobian(2 * j + 1, 2 * i,     -dFi_dxi[1][0]);
      solver_.addElementToJacobian(2 * j + 1, 2 * i + 1, -dFi_dxi[1][1]);

      solver_.addElementToJacobian(2 * j,     2 * j,     -dFi_dxj[0][0]);
      solver_.addElementToJacobian(2 * j,     2 * j + 1, -dFi_dxj[0][1]);
      solver_.addElementToJacobian(2 * j + 1, 2 * j,     -dFi_dxj[1][0]);
      solver_.addElementToJacobian(2 * j + 1, 2 * j + 1, -dFi_dxj[1][1]);
    }

 }

 if (mouse_spring_.active) {
    Particle& p0 = body_.particles[ mouse_spring_.particle_index ];

    if (!p0.locked) {
        vec2 pos0 = p0.position;
        vec2 pos1 = mouse_spring_.mouse_position;

        vec2 deltaX = pos0 - pos1;
        float normDX = norm(deltaX);

        float topLeft     = -spring_stiffness_ - spring_damping_ * deltaX[0] * deltaX[0] / (time_step_ * normDX * normDX);
        float bottomRight = -spring_stiffness_ - spring_damping_ * deltaX[1] * deltaX[1] / (time_step_ * normDX * normDX);
        float other = -spring_stiffness_ - spring_damping_ * deltaX[0] * deltaX[1] / (time_step_ * normDX * normDX);

        int i = p0.id;
        solver_.addElementToJacobian(2 * i,     2 * i,      topLeft);
        solver_.addElementToJacobian(2 * i,     2 * i + 1,  other);
        solver_.addElementToJacobian(2 * i + 1, 2 * i,      other);
        solver_.addElementToJacobian(2 * i + 1, 2 * i + 1,  bottomRight);
    }
}

 if (area_forces_) {
    std::vector<Triangle>& triangles = body_.triangles;
    for (std::vector<Triangle>::iterator triangle = triangles.begin(); triangle != triangles.end(); ++triangle) {

        const int i = triangle->particle0->id;
        const int j = triangle->particle1->id;
        const int k = triangle->particle2->id;

        const vec2& pi = triangle->particle0->position;
        const vec2& pj = triangle->particle1->position;
        const vec2& pk = triangle->particle2->position;

        const float Area = triangle->area();
        const vec2 dArea_dxi( (pj[1] - pk[1]) / 2.0f, (pk[0] - pj[0]) / 2.0f);
        const vec2 dArea_dxj( (pk[1] - pi[1]) / 2.0f, (pi[0] - pk[0]) / 2.0f);
        const vec2 dArea_dxk( (pi[1] - pj[1]) / 2.0f, (pj[0] - pi[0]) / 2.0f);

        const float C = area_stiffness_ * (triangle->area() - triangle->rest_area);
        const vec2 dC_dxi = area_stiffness_ * dArea_dxi;
        const vec2 dC_dxj = area_stiffness_ * dArea_dxj;
        const vec2 dC_dxk = area_stiffness_ * dArea_dxk;


        if (!(triangle->particle0->locked)) {
          // dFi/dxi
          solver_.addElementToJacobian(2 * i,     2 * i,      -dC_dxi[0] * dArea_dxi[0]);
          solver_.addElementToJacobian(2 * i,     2 * i + 1,  -dC_dxi[1] * dArea_dxi[0]);
          solver_.addElementToJacobian(2 * i + 1, 2 * i,      -dC_dxi[0] * dArea_dxi[1]);
          solver_.addElementToJacobian(2 * i + 1, 2 * i + 1,  -dC_dxi[1] * dArea_dxi[1]);

          // dFi/dxj
          solver_.addElementToJacobian(2 * i,     2 * j,      -dC_dxj[0] * dArea_dxi[0]);
          solver_.addElementToJacobian(2 * i,     2 * j + 1,  -dC_dxj[1] * dArea_dxi[0] - C / 2.0f);
          solver_.addElementToJacobian(2 * i + 1, 2 * j,      -dC_dxj[0] * dArea_dxi[1] + C / 2.0f);
          solver_.addElementToJacobian(2 * i + 1, 2 * j + 1,  -dC_dxj[1] * dArea_dxi[1]);

          // dFi/dxk
          solver_.addElementToJacobian(2 * i,     2 * k,      -dC_dxk[0] * dArea_dxi[0]);
          solver_.addElementToJacobian(2 * i,     2 * k + 1,  -dC_dxk[1] * dArea_dxi[0] + C / 2.0f);
          solver_.addElementToJacobian(2 * i + 1, 2 * k,      -dC_dxk[0] * dArea_dxi[1] - C / 2.0f);
          solver_.addElementToJacobian(2 * i + 1, 2 * k + 1,  -dC_dxk[1] * dArea_dxi[1]);
        }

        if (!(triangle->particle1->locked)) {
          // dFj/dxi
          solver_.addElementToJacobian(2 * j,     2 * i,      -dC_dxi[0] * dArea_dxj[0]);
          solver_.addElementToJacobian(2 * j,     2 * i + 1,  -dC_dxi[1] * dArea_dxj[0] + C / 2.0f);
          solver_.addElementToJacobian(2 * j + 1, 2 * i,      -dC_dxi[0] * dArea_dxj[1] - C / 2.0f);
          solver_.addElementToJacobian(2 * j + 1, 2 * i + 1,  -dC_dxi[1] * dArea_dxj[1]);

          // dFj/dxj
          solver_.addElementToJacobian(2 * j,     2 * j,      -dC_dxj[0] * dArea_dxj[0]);
          solver_.addElementToJacobian(2 * j + 1, 2 * j,      -dC_dxj[0] * dArea_dxj[1]);
          solver_.addElementToJacobian(2 * j,     2 * j + 1,  -dC_dxj[1] * dArea_dxj[0]);
          solver_.addElementToJacobian(2 * j + 1, 2 * j + 1,  -dC_dxj[1] * dArea_dxj[1]);

          // dFj/dxk
          solver_.addElementToJacobian(2 * j,     2 * k,      -dC_dxk[0] * dArea_dxj[0]);
          solver_.addElementToJacobian(2 * j,     2 * k + 1,  -dC_dxk[1] * dArea_dxj[0] - C / 2.0f);
          solver_.addElementToJacobian(2 * j + 1, 2 * k,      -dC_dxk[0] * dArea_dxj[1] + C / 2.0f);
          solver_.addElementToJacobian(2 * j + 1, 2 * k + 1,  -dC_dxk[1] * dArea_dxj[1]);
        }

        if (!(triangle->particle2->locked)) {
          // dFk/dxi
          solver_.addElementToJacobian(2 * k,     2 * i,      -dC_dxi[0] * dArea_dxk[0]);
          solver_.addElementToJacobian(2 * k,     2 * i + 1,  -dC_dxi[1] * dArea_dxk[0] - C / 2.0f);
          solver_.addElementToJacobian(2 * k + 1, 2 * i,      -dC_dxi[0] * dArea_dxk[1] + C / 2.0f);
          solver_.addElementToJacobian(2 * k + 1, 2 * i + 1,  -dC_dxi[1] * dArea_dxk[1]);

          // dFk/dxj
          solver_.addElementToJacobian(2 * k,     2 * j,      -dC_dxj[0] * dArea_dxk[0]);
          solver_.addElementToJacobian(2 * k,     2 * j + 1,  -dC_dxj[1] * dArea_dxk[0] + C / 2.0f);
          solver_.addElementToJacobian(2 * k + 1, 2 * j,      -dC_dxj[0] * dArea_dxk[1] - C / 2.0f);
          solver_.addElementToJacobian(2 * k + 1, 2 * j + 1,  -dC_dxj[1] * dArea_dxk[1]);

          // dFk/dxk
          solver_.addElementToJacobian(2 * k,     2 * k,      -dC_dxk[0] * dArea_dxk[0]);
          solver_.addElementToJacobian(2 * k + 1, 2 * k,      -dC_dxk[0] * dArea_dxk[1]);
          solver_.addElementToJacobian(2 * k,     2 * k + 1,  -dC_dxk[1] * dArea_dxk[0]);
          solver_.addElementToJacobian(2 * k + 1, 2 * k + 1,  -dC_dxk[1] * dArea_dxk[1]);
        }

    }
 }


}

//=============================================================================
void ImplicitSolver::solve (float dt, float mass,
                            std::vector<Particle> &particles)
{
  int num_particles = particles.size ();

  /// Build the Jacobian matrix from the kesparse set of elements
  Eigen::SparseMatrix<float> J (2 * num_particles, 2 * num_particles);
  J.setFromTriplets (triplets_.begin (), triplets_.end ());

  /// Build up the position, velocity and force vectors
  Eigen::VectorXf pos_vec (Eigen::VectorXf::Zero (2 * num_particles)),
                  velocity_vec (Eigen::VectorXf::Zero (2 * num_particles)),
                  force_vec (Eigen::VectorXf::Zero (2 * num_particles));
  for (size_t p_i = 0; p_i < num_particles; ++p_i)
  {
    pos_vec (2 * p_i + 0) = particles[p_i].position[0];
    pos_vec (2 * p_i + 1) = particles[p_i].position[1];
    velocity_vec (2 * p_i + 0) = particles[p_i].velocity[0];
    velocity_vec (2 * p_i + 1) = particles[p_i].velocity[1];
    force_vec (2 * p_i + 0) = particles[p_i].force[0];
    force_vec (2 * p_i + 1) = particles[p_i].force[1];
  }

  /// Kick out the fixed particles by creating a sparse selection matrix
  std::vector<Eigen::Triplet<float> > triplets_selection;
  int valid_particle_index = 0;
  for (size_t p_i = 0; p_i < num_particles; ++p_i)
    if (!particles[p_i].locked)
    {
      triplets_selection.push_back (Eigen::Triplet<float> (2 * p_i + 0, 2 * valid_particle_index + 0, 1.f));
      triplets_selection.push_back (Eigen::Triplet<float> (2 * p_i + 1, 2 * valid_particle_index + 1, 1.f));
      valid_particle_index ++;
    }
  Eigen::SparseMatrix<float> mat_selection (2 * num_particles, 2 * valid_particle_index);
  mat_selection.setFromTriplets (triplets_selection.begin (), triplets_selection.end ());

  /// Sparse identity matrix
  Eigen::SparseMatrix<float> Id (2 * valid_particle_index, 2 * valid_particle_index);
  Id.setIdentity ();

  /// Apply the selection matrix on each vector and the Jacobian matrix
  pos_vec = mat_selection.transpose () * pos_vec;
  velocity_vec = mat_selection.transpose () * velocity_vec;
  force_vec = mat_selection.transpose () * force_vec;
  J = mat_selection.transpose () * J * mat_selection;

  /// Build the right and left hand sides of the linear system
  Eigen::SparseMatrix<float> A = Id - dt * dt / mass * J;
  Eigen::VectorXf b;
  b = dt * velocity_vec + dt * dt / mass * force_vec + (Id - dt * dt / mass * J) * pos_vec;

  /// Solve the system and use the selection matrix again to arrange the new positions in a vector
  linear_solver_.analyzePattern (A);
  linear_solver_.compute (A);
  Eigen::VectorXf new_pos = mat_selection * linear_solver_.solve (b);

  /// Extract the positions from the solution vector and set the new positions and velocities inside the particle structures
  for (size_t p_i = 0; p_i < num_particles; ++p_i)
  {
    if (!particles[p_i].locked)
    {
      vec2 pos_update (new_pos (2 * p_i + 0), new_pos (2 * p_i + 1));
      particles[p_i].velocity = (pos_update - particles[p_i].position) / dt;
      particles[p_i].position = pos_update;
    }
  }
}
