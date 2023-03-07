//
// Created by elias on 20/02/2023.
//
#define _USE_MATH_DEFINES
#include "D3LSystem.h"
using namespace img;

Matrix D3LSystem::scaleFigure(const double scale) {
    Matrix m;
    for (int i = 1; i < 4; ++i)
        m(i,i) = scale;
    return m;
}

Matrix D3LSystem::rotateX(const double angle) {
    Matrix m;
    m(2, 2) = cos(angle);
    m(2, 3) = sin(angle);
    m(3, 3) = cos(angle);
    m(3, 2) = -sin(angle);
    return m;
}

Matrix D3LSystem::rotateY(const double angle) {
    Matrix m;
    m(1, 1) = cos(angle);
    m(1, 3) = -sin(angle);
    m(3, 1) = sin(angle);
    m(3, 3) = cos(angle);
    return m;
}

Matrix D3LSystem::rotateZ(const double angle) {
    Matrix m;
    m(1, 1) = cos(angle);
    m(1, 2) = sin(angle);
    m(2, 1) = -sin(angle);
    m(2, 2) = cos(angle);
    return m;
}

Matrix D3LSystem::translate(const Vector3D &vector) {
    Matrix m;
    m(4, 1) = vector.x;
    m(4, 2) = vector.y;
    m(4, 3) = vector.z;
    return m;
}

void D3LSystem::toPolar(const Vector3D &point, double &theta, double &phi, double &r){
    r = sqrt(point.x*point.x + point.y*point.y + point.z*point.z);
    theta = atan2(point.y, point.x);
    phi = acos(point.z / r);
}

void D3LSystem::applyTransformation(Figure &fig, const Matrix &m) {
    for(Vector3D& point: fig.points)
        point *= m;
}

Matrix D3LSystem::eyePointTrans(const Vector3D &eyepoint) {
    double theta, phi, r;
    toPolar(eyepoint, theta, phi, r);
    Matrix m = rotateZ(-1*M_PI / 2 - theta) * rotateX(-phi) * translate(Vector3D::point(0, 0, -r));
    return m;
}

Point2D D3LSystem::doProjection(const Vector3D &point, const double d) {
    double x = d * point.x / -point.z;
    double y = d * point.y / -point.z;

    return Point2D(x,y,point.z);
}

Lines2D D3LSystem::doProjection() {
    Lines2D lines;

    for (auto i : figures){
        for(auto j : i.faces){
            for(int k = 1; k < j.point_indexes.size(); k++){
                Point2D p = doProjection(i.points[j.point_indexes[k - 1]], 1);
                Point2D p2 = doProjection(i.points[j.point_indexes[k]], 1);
                Line2D lijn(p, p2, i.color);
                lines.push_back(lijn);
            }
            if(j.point_indexes.size() > 2){
                Point2D p = doProjection(i.points[j.point_indexes[j.point_indexes.size() - 1]], 1);
                Point2D p1 = doProjection(i.points[j.point_indexes[0]], 1);
                Line2D lijn(p, p1, i.color);
                lines.push_back(lijn);
            }
        }
    }
    return lines;
}

void D3LSystem::createFigure(const ini::Configuration &configuration, string &name, const Matrix& eyeMatrix) {
    int pointCount = configuration[name]["nrPoints"].as_int_or_default(0);
    int lineCount = configuration[name]["nrLines"].as_int_or_default(0);
    string type = configuration[name]["type"].as_string_or_default("");
    vector<double> center = configuration[name]["center"].as_double_tuple_or_default({0, 0, 0});
    vector<double> color = configuration[name]["color"].as_double_tuple_or_default(vector<double>({1, 1, 1}));
    double centerX = center[0];
    double centerY = center[1];
    double centerZ = center[2];
    Figure figure;
    if (type == "LineDrawing"){
        int nrPoints = configuration[name]["nrPoints"].as_int_or_die();
        int nrLines = configuration[name]["nrLines"].as_int_or_die();
        for (int i = 0; i < nrPoints; ++i) {
            string pName = "point" + to_string(i);
            vector<double> v = configuration[name][pName].as_double_tuple_or_die();
            figure.points.emplace_back(Vector3D::point(v[0], v[1], v[2]));
        }

        for (int j = 0; j < nrLines; ++j) {
            string lName = "line" + to_string(j);
            vector<int> v = configuration[name][lName].as_int_tuple_or_die();
            figure.faces.emplace_back(Face(v));
        }
    }
    else if (type == "Cube") {
        figure = createCube();
    } else if (type == "Tetrahedron") {
        figure = createTetrahedron();
    } else if (type == "Octahedron") {
        figure = createOctahedron();
    } else if (type == "Icosahedron") {
        figure = createIcosahedron();
    } else if (type == "Dodecahedron") {
        figure = createDodecahedron();
    } else if (type == "Sphere") {
        int n = configuration[name]["n"].as_int_or_die();
        figure = createSphere(0, n);
    } else if (type == "Cone") {
        int n = configuration[name]["n"].as_int_or_die();
        double height = configuration[name]["height"].as_double_or_die();
        figure = createCone(n, height);
    } else if (type == "Cylinder") {
        int n = configuration[name]["n"].as_int_or_die();
        double height = configuration[name]["height"].as_double_or_die();
        figure = createCylinder(n, height);
    } else if (type == "Torus") {
        int n = configuration[name]["n"].as_int_or_die();
        int m = configuration[name]["m"].as_int_or_die();
        double r = configuration[name]["r"].as_double_or_die();
        double R = configuration[name]["R"].as_double_or_die();
        figure = createTorus(r, R, n, m);
    }
    else if (type == "3DLSystem") {
        string fileName = configuration[name]["inputfile"].as_string_or_die();
        LParser::LSystem3D parser;
        ifstream input_stream(fileName);
        input_stream >> parser;
        input_stream.close();
        figure = create3DLsystem(parser, color);
    }
    for (int i = 0; i < pointCount; i++) {
        vector<double> point = configuration[name]["point" + to_string(i)].as_double_tuple_or_die();
        double x = point[0];
        double y = point[1];
        double z = point[2];
        figure.points.push_back(Vector3D::point(x, y, z));
    }

    for (int i = 0; i < lineCount; i++) {
        vector<int> line = configuration[name]["line" + to_string(i)].as_int_tuple_or_die();
        Face face;
        for (int j: line) {
            face.point_indexes.push_back(j);
        }
        figure.faces.push_back(face);
    }
    double xRotate = (configuration[name]["rotateX"].as_double_or_default(0)) * M_PI / 180;
    double yRotate = (configuration[name]["rotateY"].as_double_or_default(0)) * M_PI / 180;
    double zRotate = (configuration[name]["rotateZ"].as_double_or_default(0)) * M_PI / 180;
    double scale = configuration[name]["scale"].as_double_or_default(1.0);
    Matrix transf = scaleFigure(scale) * rotateX(xRotate) * rotateY(yRotate) * rotateZ(zRotate) * translate(Vector3D::point(centerX, centerY, centerZ));
    figure.color = Color(color[0], color[1], color[2]);
    applyTransformation(figure, transf);
    applyTransformation(figure, eyeMatrix);
    figures.push_back(figure);
}
string D3LSystem::getReplacements(const string &str, const LParser::LSystem3D &parser, int iterations) {
    if (iterations == 0)
        return str;

    string newString;

    for (char j: str) {
        if (j == '+' or j == '-' or j == '(' or j == ')' or j == '^' or j == '&' or j == '/' or j == '\\' or j == '|')
            newString += j;
        else
            newString += parser.get_replacement(j);
    }
    iterations--;
    return getReplacements(newString, parser, iterations);
}

Figure D3LSystem::create3DLsystem(const LParser::LSystem3D &parser, const vector<double>& color) {
    string strings = getReplacements(parser.get_initiator(), parser, parser.get_nr_iterations());
    double angle = parser.get_angle() * M_PI / 180;
    Figure figure;
    figure.color = Color(color[0], color[1], color[2]);
    Vector3D H = Vector3D::vector(1,0,0);
    Vector3D L = Vector3D::vector(0, 1, 0);
    Vector3D U = Vector3D::vector(0, 0, 1);
    stack<vector<Vector3D>> stack;
    Vector3D pos = Vector3D::point(0,0,0);
    for(char j: strings){
        if(j == '+'){
            Vector3D newH = H * cos(angle) + L * sin(angle);
            Vector3D newL = -H * sin(angle) + L * cos(angle);
            H = newH;
            L = newL;
        } else if(j == '-'){
            Vector3D newH = H * cos(-angle) + L * sin(-angle);
            Vector3D newL = -H * sin(-angle) + L * cos(-angle);
            H = newH;
            L = newL;
        } else if (j == '^') {
            Vector3D newH = H * cos(angle) + U * sin(angle);
            Vector3D newU = -H * sin(angle) + U * cos(angle);
            H = newH;
            U = newU;
        } else if (j == '&') {
            Vector3D newH = H * cos(-angle) + U * sin(-angle);
            Vector3D newU = -H * sin(-angle) + U * cos(-angle);
            H = newH;
            U = newU;
        } else if (j == '\\') {
            Vector3D newL = L * cos(angle) - U * sin(angle);
            Vector3D newU = L * sin(angle) + U * cos(angle);
            L = newL;
            U = newU;
        } else if (j == '/') {
            Vector3D newL = L * cos(-angle) - U * sin(-angle);
            Vector3D newU = L * sin(-angle) + U * cos(-angle);
            L = newL;
            U = newU;
        } else if (j == '|') {
            H = -H;
            L = -L;
        } else if(j=='(') {
            stack.push({pos, H, L, U});
        } else if(j == ')') {
            vector<Vector3D> prevPos = stack.top();
            pos = prevPos[0];
            H = prevPos[1];
            L = prevPos[2];
            U = prevPos[3];
            stack.pop();
        } else {
            Vector3D prevPos = pos;
            pos += H;
            if (parser.draw(j)) {
                figure.points.push_back(prevPos);
                figure.points.push_back(pos);
                int point1 = figure.points.size() - 1;
                int point2 = figure.points.size() - 2;
                figure.faces.emplace_back(Face({point1, point2}));
            }
        }
    }
    return figure;
}

Figure D3LSystem::createCube() {
    Figure figure;
    figure.points = {Vector3D::point(1,-1,-1),Vector3D::point(-1,1,-1),Vector3D::point(1,1,1),Vector3D::point(-1,-1,1),Vector3D::point(1,1,-1),Vector3D::point(-1,-1,-1),Vector3D::point(1,-1,1),Vector3D::point(-1,1,1)};
    figure.faces.emplace_back(Face({0, 4, 2, 6}));
    figure.faces.emplace_back(Face({4, 1, 7, 2}));
    figure.faces.emplace_back(Face({1, 5, 3, 7}));
    figure.faces.emplace_back(Face({5, 0, 6, 3}));
    figure.faces.emplace_back(Face({6, 2, 7, 3}));
    figure.faces.emplace_back(Face({0, 5, 1, 4}));
    return figure;
}

Figure D3LSystem::createTetrahedron() {
    Figure figure;
    figure.points = {Vector3D::point(1,-1,-1), Vector3D::point(-1,1,-1), Vector3D::point(1,1,1), Vector3D::point(-1,-1,1)};
    figure.faces.emplace_back(Face({0, 1, 2}));
    figure.faces.emplace_back(Face({1, 3, 2}));
    figure.faces.emplace_back(Face({0, 3, 1}));
    figure.faces.emplace_back(Face({0, 2, 3}));
    return figure;
}

Figure D3LSystem::createOctahedron() {
    Figure figure;
    figure.points = {Vector3D::point(1,0,0),Vector3D::point(0,1,0),Vector3D::point(-1,0,0),Vector3D::point(0,-1,0),Vector3D::point(0,0,-1),Vector3D::point(0,0,1)};
    figure.faces.emplace_back(Face({0, 1, 5}));
    figure.faces.emplace_back(Face({1, 2, 5}));
    figure.faces.emplace_back(Face({2, 3, 5}));
    figure.faces.emplace_back(Face({3, 0, 5}));
    figure.faces.emplace_back(Face({1, 0, 4}));
    figure.faces.emplace_back(Face({2, 1, 4}));
    figure.faces.emplace_back(Face({3, 2, 4}));
    figure.faces.emplace_back(Face({0, 3, 4}));
    return figure;
}

Figure D3LSystem::createIcosahedron() {
    Figure figure;
    figure.points.emplace_back(Vector3D::point(0,0,(sqrt(5)/2)));
    for(int i = 2; i <=6; i++){
        figure.points.emplace_back(Vector3D::point(cos(((i-2) * 2 * M_PI/5)),sin((i - 2) * 2 * M_PI / 5),0.5));
    }
    for(int i = 7; i <=11; i++){
        figure.points.emplace_back(Vector3D::point(cos(M_PI/5 + (i - 7) * 2 * M_PI/5),sin(M_PI/5 + (i - 7) * 2 * M_PI / 5), -0.5));
    }
    figure.points.emplace_back(Vector3D::point(0,0,(-sqrt(5)/2)));
    figure.faces.emplace_back(vector<int>({0,1,2}));
    figure.faces.emplace_back(vector<int>({0,2,3}));
    figure.faces.emplace_back(vector<int>({0,3,4}));
    figure.faces.emplace_back(vector<int>({0,4,5}));
    figure.faces.emplace_back(vector<int>({0,5,1}));
    figure.faces.emplace_back(vector<int>({1,6,2}));
    figure.faces.emplace_back(vector<int>({2,6,7}));
    figure.faces.emplace_back(vector<int>({2,7,3}));
    figure.faces.emplace_back(vector<int>({3,7,8}));
    figure.faces.emplace_back(vector<int>({3,8,4}));
    figure.faces.emplace_back(vector<int>({4,8,9}));
    figure.faces.emplace_back(vector<int>({4,9,5}));
    figure.faces.emplace_back(vector<int>({5,9,10}));
    figure.faces.emplace_back(vector<int>({5,10,1}));
    figure.faces.emplace_back(vector<int>({1,10,6}));
    figure.faces.emplace_back(vector<int>({11,7,6}));
    figure.faces.emplace_back(vector<int>({11,8,7}));
    figure.faces.emplace_back(vector<int>({11,9,8}));
    figure.faces.emplace_back(vector<int>({11,10,9}));
    figure.faces.emplace_back(vector<int>({11,6,10}));
    return figure;
}

Figure D3LSystem::createDodecahedron() {
    Figure figure;
    Figure icosahedron = createIcosahedron();
    for (auto face: icosahedron.faces) {
        double newX = (icosahedron.points[face.point_indexes[0]].x + icosahedron.points[face.point_indexes[1]].x + icosahedron.points[face.point_indexes[2]].x) / 3;
        double newY = (icosahedron.points[face.point_indexes[0]].y + icosahedron.points[face.point_indexes[1]].y + icosahedron.points[face.point_indexes[2]].y) / 3;
        double newZ = (icosahedron.points[face.point_indexes[0]].z + icosahedron.points[face.point_indexes[1]].z + icosahedron.points[face.point_indexes[2]].z) / 3;
        figure.points.emplace_back(Vector3D::point(newX, newY, newZ));
    }
    figure.faces.emplace_back(Face({0, 1, 2, 3, 4}));
    figure.faces.emplace_back(Face({0, 5, 6, 7, 1}));
    figure.faces.emplace_back(Face({1, 7, 8, 9, 2}));
    figure.faces.emplace_back(Face({2, 9, 10, 11, 3}));
    figure.faces.emplace_back(Face({3, 11, 12, 13, 4}));
    figure.faces.emplace_back(Face({4, 13, 14, 5, 0}));
    figure.faces.emplace_back(Face({19, 18, 17, 16, 15}));
    figure.faces.emplace_back(Face({19, 14, 13, 12, 18}));
    figure.faces.emplace_back(Face({18, 12, 11, 10, 17}));
    figure.faces.emplace_back(Face({17, 10, 9, 8, 16}));
    figure.faces.emplace_back(Face({16, 8, 7, 6, 15}));
    figure.faces.emplace_back(Face({15, 6, 5, 14, 19}));
    return figure;
}

Figure D3LSystem::createCylinder(int n, const double h) {
    Figure figure;
    for (int i = 0; i < n; ++i) {
        figure.points.push_back(Vector3D::point(cos(2 * i * M_PI / n), sin(2 * i * M_PI / n), h));
        figure.points.push_back(Vector3D::point(cos(2 * i * M_PI / n), sin(2 * i * M_PI / n), 0));
    }
    for (int j = 0; j < 2 * n; j += 2) {
        figure.faces.push_back(Face({j % (2 * n), (j + 1) % (2 * n), (j + 3) % (2 * n), (j + 2) % (2 * n)}));
    }

    vector<int> bottom;
    vector<int> top;
    for (int k = 0; k < 2 * n; k += 2) {
        top.push_back(k);
        bottom.push_back(2 * n - (k + 1));
    }
    figure.faces.emplace_back(bottom);
    figure.faces.emplace_back(top);
    return figure;
}

Figure D3LSystem::createCone(const int n, const double h) {
    Figure figure;
    for (int i = 0; i < n; ++i) {
        figure.points.push_back(Vector3D::point(cos(2 * i * M_PI / n), sin(2 * i * M_PI / n), 0));
    }
    figure.points.push_back(Vector3D::point(0, 0, h));
    for (int j = 0; j < n; ++j) {
        figure.faces.push_back(Face({n, j % n, (j + 1) % n}));
    }
    std::vector<int> bottom;
    for (int k = 0; k < n; ++k) {
        bottom.push_back(n - k - 1);
    }
    figure.faces.emplace_back(bottom);

    return figure;
}

Vector3D D3LSystem::rescalePoints(Vector3D &point){
    double r = sqrt(pow(point.x, 2) + pow(point.y, 2) + pow(point.z, 2));
    point = Vector3D::point(point.x/r,point.y/r,point.z/r);
    return point;
}

Figure D3LSystem::createSphere(const double radius, const int n) {
    Figure figure;
    Figure icosahedron = createIcosahedron();
    if(n == 0)
        figure = createIcosahedron();
    for(int i = 0; i < n; i++) {
        figure.faces.clear();
        figure.points.clear();
        for (auto it = icosahedron.faces.begin(); it != icosahedron.faces.end(); it++) {

            figure.points.push_back(icosahedron.points[(*it).point_indexes[0]]);
            figure.points.push_back(icosahedron.points[(*it).point_indexes[1]]);
            figure.points.push_back(icosahedron.points[(*it).point_indexes[2]]);

            figure.points.push_back(Vector3D::point(((icosahedron.points[(*it).point_indexes[0]].x + icosahedron.points[(*it).point_indexes[1]].x) / 2),
                                                     ((icosahedron.points[(*it).point_indexes[0]].y + icosahedron.points[(*it).point_indexes[1]].y) / 2),
                                                     ((icosahedron.points[(*it).point_indexes[0]].z + icosahedron.points[(*it).point_indexes[1]].z) / 2)));
            figure.points.push_back(Vector3D::point(((icosahedron.points[(*it).point_indexes[0]].x + icosahedron.points[(*it).point_indexes[2]].x) / 2),
                                                     ((icosahedron.points[(*it).point_indexes[0]].y + icosahedron.points[(*it).point_indexes[2]].y) / 2),
                                                     ((icosahedron.points[(*it).point_indexes[0]].z + icosahedron.points[(*it).point_indexes[2]].z) / 2)));
            figure.points.push_back(Vector3D::point(((icosahedron.points[(*it).point_indexes[1]].x + icosahedron.points[(*it).point_indexes[2]].x) / 2),
                                                     ((icosahedron.points[(*it).point_indexes[1]].y + icosahedron.points[(*it).point_indexes[2]].y) / 2),
                                                     ((icosahedron.points[(*it).point_indexes[1]].z + icosahedron.points[(*it).point_indexes[2]].z) / 2)));
            int size = figure.points.size();
            int A = size - 6;
            int B = size - 5;
            int C = size - 4;
            int D = size - 3;
            int E = size - 2;
            int F = size - 1;

            figure.faces.emplace_back(Face({A, D, E}));
            figure.faces.emplace_back(Face({B, F, D}));
            figure.faces.emplace_back(Face({C, E, F}));
            figure.faces.emplace_back(Face({D, F, E}));
        }
        icosahedron.points = figure.points;
        icosahedron.faces = figure.faces;
        icosahedron.color = figure.color;
    }

    for(auto it = figure.points.begin(); it < figure.points.end(); it++){
        rescalePoints(*it);
    }
    return figure;
}
int findIndex(int i, int j, int m) {
    return i * m + j;
}
Figure D3LSystem::createTorus(double r, double R, int n, int m) {
    Figure figure;
    for (int i = 0; i < n; ++i) {
        double u = 2 * i * M_PI / n;
        for (int j = 0; j < m; j++) {
            double v = 2 * j * M_PI / m;
            figure.points.emplace_back(Vector3D::point((R + r * cos(v)) * cos(u),
                                                      (R + r * cos(v)) * sin(u),
                                                      r * sin(v)));
            figure.faces.emplace_back(Face({i*m + j, (((i + 1)%n )* m) + j, (((i + 1)%n )* m) + (j + 1)%m, (i * m)  + (j + 1)%m}));
        }
    }
    return figure;
}

Faces3D D3LSystem::triangulate(const Face &face) {
    Faces3D triangles;
    int beginPoint = face.point_indexes[0];
    for (int i = 1; i < face.point_indexes.size() - 1; ++i) {
        triangles.emplace_back(vector<int>({beginPoint, face.point_indexes[i], face.point_indexes[i+1]}));
    }
    return triangles;
}

