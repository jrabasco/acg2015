TARGET = acg-physics
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

#QT += opengl

INCLUDEPATH += src

HEADERS += src/Mass_spring_system.h \
           src/Mass_spring_viewer.h \
           src/utils/gl.h \
           src/utils/gl_text.h \
           src/utils/GLUT_viewer.h \
           src/utils/vec2.h \
           src/utils/Viewer_2D.h

SOURCES += src/main.cpp \
           src/Mass_spring_system.cpp \
           src/Mass_spring_viewer.cpp \
           src/utils/gl_text.cpp \
           src/utils/GLUT_viewer.cpp \
           src/utils/Viewer_2D.cpp

INCLUDEPATH += $$PWD/external

unix:macx {
    QMAKE_CXXFLAGS_X86_64 += -mmacosx-version-min=10.9 -stdlib=libc++
    QMAKE_LFLAGS += -framework GLUT
    QMAKE_LFLAGS += -framework opengl
    QMAKE_CXXFLAGS += -Wno-deprecated-declarations -stdlib=libc++
}
unix:!macx {
    LIBS += -lGL -lGLU -lglut
}

QMAKE_CXXFLAGS += -Wall -Wno-unused-parameter
QMAKE_CXXFLAGS_RELEASE += -DNDEBUG
