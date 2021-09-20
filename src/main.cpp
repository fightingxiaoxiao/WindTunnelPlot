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

/*
-----------输入---------------
参数 1: 部件几何 (.ply)
参数 2: 部件原始数据点 (.csv)

-----------输出---------------
参数 3: 插值后的部件数据文件(.csv)
参数 4: 部件网格文件(.csv)

-----------其他---------------
参数 5: 迭代次数
参数 6: block在x方向的剖分
参数 7: block在y方向的剖分
*/

#include "mesh.h"

int main(int argc, char *argv[])
{
    Word   filename = (Word)argv[1];
    Word   data_path = (Word)argv[2];
    Word   outfile_value = (Word)argv[3];
    Word   outfile_cell = (Word)argv[4];
    Word   zone_avg_file = (Word)argv[5];
    size_t iterations = std::stoi(argv[6]);
    size_t block_x = std::stoi(argv[7]);
    size_t block_y = std::stoi(argv[8]);
    // Word filename = "F:/wind_tunnel_data/data/AU.ply";
    // Word data_path = "F:/wind_tunnel_data/data/AU.csv";
    SMILE::Mesh       mesh;
    SMILE::pointCloud point_cloud;
    /*
       if (mesh.readPLY(filename))
       {
           std::cout << "Cannot find file: " << filename << std::endl;
           return EXIT_FAILURE;
       }
    */
    if (mesh.readOBJGroup(filename))
    {
        std::cout << "failed to read file: " << filename << std::endl;
        return EXIT_FAILURE;
    }
    mesh.getMaxMin();

    if (point_cloud.readData(data_path, 0))
    {
        std::cout << "Cannot find file: " << data_path << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "xmax: " << mesh.xmax << std::endl;
    std::cout << "xmin: " << mesh.xmin << std::endl;
    std::cout << "ymax: " << mesh.ymax << std::endl;
    std::cout << "ymin: " << mesh.ymin << std::endl;

    Scalar           x_offset = (mesh.xmax - mesh.xmin) * 0.05;
    Scalar           y_offset = (mesh.ymax - mesh.ymin) * 0.05;
    SMILE::blockMesh block = SMILE::blockMesh(block_x, block_y, mesh.xmax + x_offset, mesh.xmin - x_offset, mesh.ymax + y_offset, mesh.ymin - y_offset);

    // std::cout << block.storage << std::endl;
    std::cout << std::endl;
    block.mapValueFromPointCloud(&point_cloud, (Word) "p");
    // std::cout << block.storage << std::endl;

    block.iter(iterations);
    std::cout << std::endl;

    mesh.mapValueFromBlockMesh(&block, (Word) "p");

    mesh.triangulation();

    mesh.interpolateScalarValueFromNodeToCell();

    // mesh.getNewCellValue((Word) "p", (Word) "p_coeff", 1 / 0.5 / 1.225 / 6.2 / 6.2);
    mesh.writeZoneAverage((Word) "p", zone_avg_file);

    mesh.writeNodeValue(outfile_value, (Word) "p");
    mesh.writeCellInfo(outfile_cell);

    // std::cout << mesh.cells.size() << std::endl;

    return EXIT_SUCCESS;
}