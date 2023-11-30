/**
 * @file    tree_mesh_builder.cpp
 *
 * @author  David Chocholaty <xchoch09@stud.fit.vutbr.cz>
 *
 * @brief   Parallel Marching Cubes implementation using OpenMP tasks + octree early elimination
 *
 * @date    30 November 2023
 **/

#include <iostream>
#include <math.h>
#include <limits>

#include "tree_mesh_builder.h"

constexpr float sqrt3Div2 = 0.86602540378f;
constexpr float gridSizeCutoff = 1.0f;

TreeMeshBuilder::TreeMeshBuilder(unsigned gridEdgeSize)
    : BaseMeshBuilder(gridEdgeSize, "Octree")
{

}

unsigned TreeMeshBuilder::treeTraversal(const Vec3_t<float> &pos,
                                        const ParametricScalarField &field,
                                        const float gridSize)
{
    unsigned totalTriangles = 0;
    const float edgeSize = gridSize * mGridResolution;
    const float halfGridSize = gridSize / 2.0f;

    Vec3_t<float> cubeCenter3D((pos.x + halfGridSize) * mGridResolution,
                               (pos.y + halfGridSize) * mGridResolution,
                               (pos.z + halfGridSize) * mGridResolution);

    if (evaluateFieldAt(cubeCenter3D, field) > mIsoLevel + sqrt3Div2 * edgeSize) {
        return 0;
    }

    // Based on the assignment the grid size values are only powers of two. That means, that after some amount of steps the gridSize value has to be one. So, we can test only equality instead of less or equal condition.
    if (gridSize == gridSizeCutoff) {
        totalTriangles += buildCube(pos, field);
    } else {
        for (int i = 0; i < 8; i++) {
            Vec3_t<float> newPos(pos.x + halfGridSize * sc_vertexNormPos[i].x,
                                 pos.y + halfGridSize * sc_vertexNormPos[i].y,
                                 pos.z + halfGridSize * sc_vertexNormPos[i].z);

            totalTriangles += treeTraversal(newPos, field, halfGridSize);
        }
    }

    return totalTriangles;
}

unsigned TreeMeshBuilder::marchCubes(const ParametricScalarField &field)
{
    // Suggested approach to tackle this problem is to add new method to
    // this class. This method will call itself to process the children.
    // It is also strongly suggested to first implement Octree as sequential
    // code and only when that works add OpenMP tasks to achieve parallelism.
    const float gridSize = static_cast<const float>(mGridSize);

    unsigned totalTriangles = treeTraversal(sc_vertexNormPos[0], field, gridSize);

    return totalTriangles;
}

float TreeMeshBuilder::evaluateFieldAt(const Vec3_t<float> &pos, const ParametricScalarField &field)
{
    // NOTE: This method is called from "buildCube(...)"!

    // 1. Store pointer to and number of 3D points in the field
    //    (to avoid "data()" and "size()" call in the loop).
    const Vec3_t<float> *pPoints = field.getPoints().data();
    const unsigned count = unsigned(field.getPoints().size());

    float value = std::numeric_limits<float>::max();

    // 2. Find minimum square distance from points "pos" to any point in the
    //    field.
    for(unsigned i = 0; i < count; ++i)
    {
        float distanceSquared  = (pos.x - pPoints[i].x) * (pos.x - pPoints[i].x);
        distanceSquared       += (pos.y - pPoints[i].y) * (pos.y - pPoints[i].y);
        distanceSquared       += (pos.z - pPoints[i].z) * (pos.z - pPoints[i].z);

        // Comparing squares instead of real distance to avoid unnecessary
        // "sqrt"s in the loop.
        value = std::min(value, distanceSquared);
    }

    // 3. Finally take square root of the minimal square distance to get the real distance
    return sqrt(value);
}

void TreeMeshBuilder::emitTriangle(const BaseMeshBuilder::Triangle_t &triangle)
{
    // NOTE: This method is called from "buildCube(...)"!

    // Store generated triangle into vector (array) of generated triangles.
    // The pointer to data in this array is return by "getTrianglesArray(...)" call
    // after "marchCubes(...)" call ends.
    #pragma omp critical
    mTriangles.push_back(triangle);
}
