//
// Created by elias on 28/02/2023.
//

#ifndef ENGINE_ZBUFFER_H
#include <limits>
#include <iostream>
#include <fstream>
#include "easy_image.h"
#define ENGINE_ZBUFFER_H

class ZBuffer: public std::vector<std::vector<double>>{
public:
    ZBuffer(const int width, const int height);
    vector<vector<double>> buf;
};
#endif //ENGINE_ZBUFFER_H
