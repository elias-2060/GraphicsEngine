//
// Created by elias on 20/02/2023.
//

#ifndef ENGINE_FACE_H
#include "vector"
#define ENGINE_FACE_H
using namespace std;
class Face{
public:
    vector<int> point_indexes;
    Face()= default;;
    explicit Face(const std::vector<int> &points) : point_indexes(points) {}
};
typedef vector<Face> Faces3D;

#endif //ENGINE_FACE_H
