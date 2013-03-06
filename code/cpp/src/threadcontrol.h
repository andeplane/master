#pragma once
#include <set>
#include <vector>
class System;

using namespace std;

class ThreadNode
{
public:
    int rank;
    set<int> owned_cells;
    set<int> connected_cells;
};

class ThreadControl
{
public:
    vector<ThreadNode> nodes;

    ThreadControl();
    void setup(int nodes, int cells_x, int cells_y, int cells_z, System *system_);

};
