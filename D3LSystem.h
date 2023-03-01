//
// Created by elias on 20/02/2023.
//

#ifndef ENGINE_D3LSYSTEM_H
#include "easy_image.h"
#include "vector3d.h"
#include "ini_configuration.h"
#include "l_parser.h"
#include "fstream"
#include "stack"
#include "Figure.h"
#define ENGINE_D3LSYSTEM_H

class D3LSystem{
private:
    Figures3D figures;
public:
    Matrix scaleFigure(const double scale);
    Matrix rotateX(const double angle);
    Matrix rotateY(const double angle);
    Matrix rotateZ(const double angle);
    Matrix translate(const Vector3D &vector);
    void applyTransformation(Figure &fig, const Matrix &m);
    Matrix eyePointTrans(const Vector3D &eyepoint);
    void toPolar(const Vector3D &point, double &theta, double &phi, double &r);
    img::Point2D doProjection(const Vector3D &point, const double d);
    img::Lines2D doProjection();
    void createFigure(const ini::Configuration &configuration, string& name, const Matrix& eyeMatrix);
    Figure create3DLsystem(const LParser::LSystem3D &parser, const vector<double>& color);
    string getReplacements(const string &str, const LParser::LSystem3D &parser, int iterations);
    Figure createCube();
    Figure createTetrahedron();
    Figure createOctahedron();
    Figure createIcosahedron();
    Figure createDodecahedron();
    Figure createCylinder(int n, const double h);
    Figure createCone(const int n, const double h);
    Figure createSphere(const double radius, const int n);
    Figure createTorus(double r, double R, int n, int m);
    Vector3D rescalePoints(Vector3D &point);
};

#endif //ENGINE_D3LSYSTEM_H
