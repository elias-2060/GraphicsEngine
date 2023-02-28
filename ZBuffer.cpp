//
// Created by elias on 28/02/2023.
//
#include "ZBuffer.h"

ZBuffer::ZBuffer(const int width, const int height) {
    buf = vector<std::vector<double>>(width, vector<double>(height, std::numeric_limits<double>::infinity()));
}
