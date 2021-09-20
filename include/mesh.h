/*
 *   Copyright (c) 2020 Chen Xiaoxiao
 *   All rights reserved.

 *   Licensed under the Apache License, Version 2.0 (the "License");
 *   you may not use this file except in compliance with the License.
 *   You may obtain a copy of the License at

 *   http://www.apache.org/licenses/LICENSE-2.0

 *   Unless required by applicable law or agreed to in writing, software
 *   distributed under the License is distributed on an "AS IS" BASIS,
 *   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *   See the License for the specific language governing permissions and
 *   limitations under the License.
 */

#ifndef MESH_H
#define MESH_H

#include "string_split.h"
#include "type.h"
#include <fstream>
#include <iostream>
#include <omp.h>
#include <cfloat>
// using namespace std;
namespace SMILE
{
class blockMesh;
class Point
{
public:
    Label  id;
    Word   name;
    Vector coord;
    Vector normal;
    Word   zone_name;

    Dict<Word, Scalar> scalar = {};
};

class pointCloud
{
public:
    List<Point> storage;

    bool readData(Word filename, size_t field_index);
};

class Node
{
public:
    Label              id;
    Vector             coord;
    Dict<Word, Scalar> scalar = {};
    Dict<Word, Vector> vector = {};
    List<Label>        cells;

    //Node();
};

class Cell
{
public:
    Label              id;
    Vector             center;
    List<Label>        nodes;
    size_t             size;
    Word               zone_name;
    Dict<Word, Scalar> scalar = {};
    Dict<Word, Vector> vector = {};
};

class Mesh
{
public:
    List<Node> nodes;
    List<Cell> cells;

    Scalar xmax = -1 * DBL_MAX;
    Scalar xmin = DBL_MAX;
    Scalar ymax = -1 * DBL_MAX;
    Scalar ymin = DBL_MAX;

    bool readPLY(Word filename);
    bool readOBJGroup(Word filename);
    bool writeVTK(Word filename, Word valueName);
    bool writeNodeValue(Word filename, Word valueName);
    bool writeCellInfo(Word filename);
    void getMaxMin();

    Dict<Word, Scalar> writeZoneAverage(Word valueType, Word filename);

    void mapValueFromPointCloud(pointCloud *, Word);
    void mapValueFromBlockMesh(blockMesh *, Word);
    void interpolateScalarValueFromNodeToCell();
    void triangulation();

    void getNewCellValue(Word existValueType, Word newValueType, Scalar coeff);
};

// 2D Octree
class Octree
{
public:
    Word          name;
    Octree *      sub_nodes[4];
    Vector        Max;
    Vector        Min;
    List<Point *> storage;
};

class blockMesh
{
public:
    Word name;

    Scalar xmax;
    Scalar xmin;
    Scalar ymax;
    Scalar ymin;

    Scalar x_delta;
    Scalar y_delta;
    size_t X;
    size_t Y;

    blockMesh(size_t, size_t, Scalar x_max, Scalar x_min, Scalar y_max, Scalar y_min);

    Eigen::MatrixXd storage;
    Eigen::MatrixXd storage_old;
    Eigen::MatrixXi storage_fix;
    void            mapValueFromPointCloud(pointCloud *, Word);
    void            iter(size_t);
};
} // namespace SMILE

#endif