#ifndef PROJECT4_UTILS_H
#define PROJECT4_UTILS_H
#include <iostream>
#include <unordered_map>
#include <string>
#include <stack>
#include <limits>
#include <vector>
using namespace std;
class Vertex {
public:
    //string name;
    int ID;
    //unordered_map<string, pair<int,int>> adjacencyMap;
    vector<pair<int,int>> adjacencyMap;
    Vertex(int ID=-1);
    void addEdge(string name, pair<int,int> capacity_cost);
    friend ostream& operator<<(ostream&output,const Vertex v);
    //ostream& operator<<(ostream& os, const Vertex v);
};

class AugmentingPath {
public:
    int amount;
    int totalCost;
    stack <string> path;
    void printPath(stack<string> path);
    AugmentingPath();
    AugmentingPath(int amount,int totalCost, stack<string> path);
};

#endif //PROJECT4_UTILS_H
