//
// Created by elias on 20/02/2023.
//

#ifndef ENGINE_FIGURE_H
#include "vector3d.h"
#include "Face.h"
#include "easy_image.h"
#define ENGINE_FIGURE_H
using namespace std;

class Figure{
public:
    vector<Vector3D> points;
    Faces3D faces;
    img::Color color;
};
typedef list<Figure> Figures3D;
#endif //ENGINE_FIGURE_H
