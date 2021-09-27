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

#include "mesh.h"
#include "csvparser.h"
namespace SMILE
{
    bool pointCloud::readData(Word filename, size_t field_index)
    {
        std::cout << "Read .csv file from " << filename << "." << std::endl;
        std::fstream _file;
        _file.open(filename.c_str(), std::ios::in);
        if (!_file)
            return EXIT_FAILURE;

        CsvParser *csvparser = CsvParser_new(filename.c_str(), ",", 0);
        CsvRow *row;

        while ((row = CsvParser_getRow(csvparser)))
        {
            const char **rowFields = CsvParser_getFields(row);

            Point p_temp;
            p_temp.name = rowFields[0];
            p_temp.coord << std::stod(rowFields[1]), std::stod(rowFields[2]), std::stod(rowFields[3]);
            p_temp.scalar[(Word) "p"] = std::stod(rowFields[field_index + 4]);

            storage.push_back(p_temp);
            CsvParser_destroy_row(row);
        }
        CsvParser_destroy(csvparser);

        return EXIT_SUCCESS;
    }

    bool Mesh::readOBJGroup(Word filename)
    {
        std::cout << "Read .obj file from " << filename << "." << std::endl;

        Word filePath = filename;
        std::ifstream file;
        file.open(filePath, std::ios::in);
        if (!file)
        {
            std::cerr << "Failed to read .obj file." << std::endl;
            return EXIT_FAILURE;
        }
        Word strLine;
        Word zone_name = (Word) "";
        // 逐行读取
        while (getline(file, strLine))
        {

            std::istringstream out(strLine);
            Word split_str;
            Node node_temp;
            Cell cell_temp;

            if (strLine.size() < 1)
                continue;

            auto line_list = splitStr(strLine, (Word) " ");

            if (line_list[0] == (Word) "o")
            {
                if (line_list.size() > 1)
                {
                    zone_name = line_list[1];
                }
                else
                {
                    std::cerr << "Cannot get object name in .obj file." << std::endl;
                    return EXIT_FAILURE;
                }
            }

            if (line_list[0] == (Word) "v")
                if (line_list.size() > 3)
                {
                    Scalar x, y, z;
                    std::stringstream x_str(line_list[1]), y_str(line_list[2]), z_str(line_list[3]);
                    x_str >> x;
                    y_str >> y;
                    z_str >> z;

                    node_temp.coord << x, y, z;
                    node_temp.id = nodes.size();
                    nodes.push_back(node_temp);
                    node_temp = Node();
                }
                else
                {
                    std::cerr << "illegal vertice info in .obj file." << std::endl;
                    return EXIT_FAILURE;
                }

            if (line_list[0] == (Word) "f")
            {
                size_t index;
                for (int i = 1; i < line_list.size(); i++)
                {
                    auto temp = splitStr(line_list[i], (Word) "/");

                    std::stringstream sin(temp[0]);
                    sin >> index;
                    if (index - 1 >= nodes.size())
                    {
                        std::cerr << "exceed node number: " << index << std::endl;
                        return EXIT_FAILURE;
                    }
                    cell_temp.nodes.push_back(index - 1);
                }
                cell_temp.zone_name = zone_name;
                cells.push_back(cell_temp);
                cell_temp = Cell();
            }
        }
        return EXIT_SUCCESS;
    }

    bool Mesh::readPLY(Word filename)
    {
        std::cout << "Read .ply file from " << filename << "." << std::endl;

        Word filePath = filename;
        std::ifstream file;
        file.open(filePath, std::ios::in);
        if (!file)
            return EXIT_FAILURE;

        float data;
        int index;
        bool not_head = false;
        Word strLine;
        Word header = "end_header";
        size_t count = 0;

        Label node_count = 0;
        Label cell_count = 0;
        // 逐行读取
        while (getline(file, strLine))
        {
            if (strLine.empty())
                continue;
            // std::cout << strLine << std::endl;
            if (strLine == header)
            {
                not_head = true;
                continue;
            }
            if (not_head)
            {
                std::istringstream out(strLine);
                Word split_str;
                bool isCellData = false;
                int item_count = 0;

                Node node_temp;
                Cell cell_temp;

                while (out >> split_str)
                {
                    std::stringstream sin(split_str);
                    if (item_count == 0 && split_str.size() == 1)
                    {
                        sin >> index;
                        isCellData = true;
                        cell_temp.size = index;
                        item_count++;
                        // std::cout << index << std::endl;
                        continue;
                    }
                    if (isCellData)
                    {
                        sin >> index;
                        cell_temp.nodes.push_back(index);
                        // std::cout << data << std::endl;

                        if (item_count == cell_temp.size)
                        {
                            cell_temp.id = cell_count;
                            cells.push_back(cell_temp);
                            cell_count++;
                        }
                    }
                    else
                    {
                        if (item_count < 3)
                        {
                            sin >> data;
                            node_temp.coord(item_count) = data;
                        }
                        if (item_count == 2)
                        {
                            node_temp.id = node_count;
                            nodes.push_back(node_temp);
                            // std::cout << node_temp.coord.transpose() <<
                            // std::endl;
                            node_count++;
                        }
                        // std::cout << data << std::endl;
                    }
                    item_count++;
                }
            }
            count++;
        }

        return EXIT_SUCCESS;
    }

    bool Mesh::writeVTK(Word filename, Word valueName)
    {
        /* header */

        return EXIT_SUCCESS;
    }

    bool Mesh::writeNodeValue(Word filename, Word valueName)
    {
        std::ofstream outfile(filename.c_str(), std::ios::trunc);
        for (auto node : nodes)
            outfile << node.coord[0] << "," << node.coord[1] << "," << node.coord[2] << "," << node.scalar[valueName] << std::endl;

        return EXIT_SUCCESS;
    }

    bool Mesh::writeCellInfo(Word filename)
    {
        std::ofstream outfile(filename.c_str(), std::ios::trunc);
        for (auto cell : cells)
        {
            int count = 0;
            for (auto i : cell.nodes)
            {
                if (count < 2)
                    outfile << nodes[i].id << ",";
                else if (count == 2)
                    outfile << nodes[i].id << std::endl;
                else
                    std::cerr << "node over 3." << std::endl;
                count++;
            }
        }
        return EXIT_SUCCESS;
    }

    void Mesh::mapValueFromPointCloud(pointCloud *data_set, Word interp_type)
    {
        if (interp_type == (Word) "RBF")
            for (auto node : nodes)
                for (auto point : data_set->storage)
                    break;
        // return EXIT_SUCCESS;
    }

    void Mesh::mapValueFromBlockMesh(blockMesh *block, Word valueType)
    {
        size_t lost_count = 0;
        for (auto &node : nodes)
        {
            if (node.coord[0] < xmin || node.coord[1] < ymin || node.coord[0] > xmax || node.coord[1] > ymax)
            {
                lost_count++;
                continue;
            }
            size_t x_i = (size_t)((node.coord[0] - block->xmin) / block->x_delta);
            size_t y_i = (size_t)((node.coord[1] - block->ymin) / block->y_delta);

            node.scalar[valueType] = block->storage(x_i, y_i);
            // std::cout << node.scalar[valueType] << std::endl;
        }

        if (lost_count > 0)
            std::cerr << "Warning: " << lost_count / nodes.size() * 100 << "%% nodes out of the block.\n";
    }

    void Mesh::interpolateScalarValueFromNodeToCell()
    {
        for (auto &cell : cells)
        {
            cell.center << 0., 0., 0.;
            for (auto i : cell.nodes)
            {
                // std::cout << nodes[i].coord.transpose() << std::endl;
                cell.center += nodes[i].coord;
            }
            cell.center /= (Scalar)cell.nodes.size();
            // std::cout << cell.center.transpose() << std::endl;
            Scalar total_weight = 0.;
            Scalar current_weight;
            for (auto i : cell.nodes)
            {
                for (auto value : nodes[i].scalar)
                {
                    current_weight = 1. / (cell.center - nodes[i].coord).squaredNorm();
                    if (cell.scalar.find(value.first) == cell.scalar.end())
                        cell.scalar[value.first] = 0.;

                    cell.scalar[value.first] = cell.scalar[value.first] * total_weight + value.second * current_weight;
                    cell.scalar[value.first] /= total_weight + current_weight;
                }
                total_weight += current_weight;
            }
        }
    }

    void Mesh::getMaxMin()
    {
        for (auto n : nodes)
        {
            xmax = (xmax >= n.coord[0]) ? xmax : n.coord[0];
            xmin = (xmin <= n.coord[0]) ? xmin : n.coord[0];

            ymax = (ymax >= n.coord[1]) ? ymax : n.coord[1];
            ymin = (ymin <= n.coord[1]) ? ymin : n.coord[1];

            // std::cerr << n.coord[0] << ", " << n.coord[1] << std::endl;
        }
    }

    void Mesh::triangulation()
    {
        List<Cell> cells_new;
        // std::cerr << cells.size() << std::endl;
        for (auto c : cells)
        {
            // std::cerr << c.size << ", " << c.id << " " << c.nodes[1]->id << " "
            //          << c.nodes[2]->id << std::endl;
            size_t cell_count = 0;
            for (int i = 0; i < c.nodes.size() - 2; i++)
            {
                Cell cell_temp;
                cell_temp.nodes.push_back(c.nodes[0]);
                cell_temp.nodes.push_back(c.nodes[i + 1]);
                cell_temp.nodes.push_back(c.nodes[i + 2]);
                cell_temp.size = 3;
                cell_temp.id = cell_count;
                cell_temp.zone_name = c.zone_name;
                cells_new.push_back(cell_temp);
                cell_count++;
            }
        }
        cells = {};
        for (auto c : cells_new)
            cells.push_back(c);
        // std::cerr << cells.size() << std::endl;
    }

    void Mesh::getNewCellValue(Word existValueType, Word newValueType, Scalar coeff)
    {
        for (auto &cell : cells)
        {
            cell.scalar[newValueType] = cell.scalar[existValueType] * coeff;
        }
    }

    Dict<Word, Scalar> Mesh::writeZoneAverage(Word valueType, Word filename)
    {
        Dict<Word, Scalar> result = {};
        Dict<Word, Scalar> weight = {};
        // Dict<Word, Scalar>::iterator iter;
        for (auto cell : cells)
        {
            if (result.find(cell.zone_name) == result.end())
            {
                // std::cout << cell.zone_name << std::endl;
                result[cell.zone_name] = 0.;
                weight[cell.zone_name] = 0.;
            }

            result[cell.zone_name] += cell.scalar[valueType];
            weight[cell.zone_name] += 1.;
        }

        for (auto &value : result)
        {
            value.second /= weight[value.first];
        }

        std::cout << "write zone average file to " << filename << "...";
        std::ofstream outfile(filename.c_str(), std::ios::trunc);
        for (auto value : result)
            outfile << value.first << "," << value.second << std::endl;
        std::cout << "done." << std::endl;
        return result;
    }

    blockMesh::blockMesh(size_t x_div, size_t y_div, Scalar x_max, Scalar x_min, Scalar y_max, Scalar y_min)
    {
        x_delta = (x_max - x_min) / x_div;
        y_delta = (y_max - y_min) / y_div;

        xmax = x_max;
        xmin = x_min;
        ymax = y_max;
        ymin = y_min;

        X = x_div + 1;
        Y = y_div + 1;
        this->storage = Eigen::MatrixXd::Constant(x_div + 1, y_div + 1, 0.);
        this->storage_fix = Eigen::MatrixXi::Constant(x_div + 1, y_div + 1, 0);
    }

    void blockMesh::mapValueFromPointCloud(pointCloud *pc, Word valueType)
    {
        size_t lost_count = 0;
        size_t count = 0;
        for (auto p : pc->storage)
        {
            if (p.coord[0] < xmin || p.coord[1] < ymin || p.coord[0] > xmax || p.coord[1] > ymax)
            {
                lost_count++;
                continue;
            }
            size_t x_i = (int)((p.coord[0] - xmin) / x_delta);
            size_t y_i = (int)((p.coord[1] - ymin) / y_delta);

            this->storage(x_i, y_i) = p.scalar[valueType];
            this->storage_fix(x_i, y_i) = 1;
            count++;
        }
        //std::cout << count << std::endl;
        // std::cout << storage_fix << std::endl;

        if (lost_count > 0)
            std::cerr << "Warning: " << lost_count << " points out of the block.\n";
    }

    void blockMesh::iter(size_t MAX_TIMES)
    {
        int n = 0;
        Scalar residual;
        std::cout << "Start Iteration..." << std::endl;
        while (n < MAX_TIMES)
        {
            residual = 0;
            int i, j;
            //this->storage_old = this->storage;
#pragma omp parallel for
            for (i = 0; i < X; i++)
            {
                Scalar residual_y = 0;
                for (j = 0; j < Y; j++)
                {
                    if (!storage_fix(i, j))
                    {
                        Scalar left_value = (i > 0) ? storage(i - 1, j) : storage(i + 1, j);
                        Scalar right_value = (i < X - 1) ? storage(i + 1, j) : storage(i - 1, j);
                        Scalar up_value = (j < Y - 1) ? storage(i, j + 1) : storage(i, j - 1);
                        Scalar down_value = (j > 0) ? storage(i, j - 1) : storage(i, j + 1);
                        Scalar last_value = this->storage(i, j);
                        this->storage(i, j) = ((left_value + right_value) / x_delta / x_delta / 2 + (down_value + up_value) / y_delta / y_delta / 2) / (1 / x_delta / x_delta + 1 / y_delta / y_delta);

                        Scalar current_residual;
                        current_residual = fabs(last_value - this->storage(i, j));
                        {
                            if (residual_y < current_residual)
                                residual_y = current_residual;
                        }
                    }
                }
#pragma omp critical
                {
                    if (residual < residual_y)
                        residual = residual_y;
                }
            }
            n++;
        }
        std::cout << "The max residual after " << MAX_TIMES << " iterations is " << residual << std::endl;
    } // namespace SMILE
} // namespace SMILE