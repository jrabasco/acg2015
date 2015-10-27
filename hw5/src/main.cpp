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

#include "utils/gl.h"
#include "Mass_spring_viewer.h"

//=============================================================================


int main(int argc, char **argv)
{
    glutInit(&argc, argv);
    Mass_spring_viewer window("Mass Springs", 1024, 768);
    glutMainLoop();
}


//=============================================================================
