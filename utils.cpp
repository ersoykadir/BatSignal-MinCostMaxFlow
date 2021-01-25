#include <iostream>
#include <unordered_map>
#include <string>
#include <stack>
#include <limits>
#include "utils.h"
using namespace std;



Vertex::Vertex(string name,int ID) {
    this->name = name;
    this->ID = ID;
}

void Vertex::addEdge(string name, pair<int,int> capacity_cost) {
    adjacencyMap[name] = capacity_cost;
}


ostream& operator<<(ostream& os, const Vertex v) {
    for (auto map: v.adjacencyMap) {
        os << map.first << " " << map.second.second << " ";
    }
    os << endl;
    return os;
}

    AugmentingPath::AugmentingPath() {
        this->amount = -1;
    }

    AugmentingPath::AugmentingPath(int amount,int totalCost, stack<string> path) {//amount is flow not cost
        this->amount = amount;
        this->totalCost = totalCost;
        this->path = path;
    }


//void printGraph(unordered_map<string, Vertex>& g) {
//    for (auto v: g) {
//        cout << v.first << " " << v.second;
//    }
//}

void printPath(stack<string> path) {
    cout << "Found path: ";
    while (path.size() > 0) {
        cout << path.top() << " ";
        path.pop();
    }
    cout << endl;
}
