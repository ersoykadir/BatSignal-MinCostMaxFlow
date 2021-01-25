#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <queue>
#include <algorithm>
#include <unordered_map>
//#include "utils.h"
#include <string>
#include <stack>
using namespace std;
//AugmentingPath findAugmentingPath(unordered_map<string, Vertex>& residualGraph, string sourceName, string targetName, string notReached);
//void augment(unordered_map<string, Vertex>& residualGraph, AugmentingPath augmentingPath);
// void printGraph(unordered_map<string, Vertex>& g) {
//     for (auto v: g) {
//         cout << v.first << " " << v.second;
//     }
// }

int sum(vector<int> *cycle,vector<unordered_map<int,int>> & residualGraph){
    //auto iter = (*cycle).begin();

    int firstV = -1;
    //string source = (*cycle)[0];
    int secondV = -1;
    int sum = 0;

    for (int i = (*cycle).size()-1; i >= 0; i--) {
        //(*sCycle).push((*cycle)[i]);
        secondV = (*cycle)[i];
        if (firstV != -1){
            int w = residualGraph[firstV][secondV];
            sum+=w;
            //secondV to firstV negative edge
            residualGraph[secondV][firstV] = w*(-1);
            residualGraph[firstV].erase(secondV);
        }
        firstV = (*cycle)[i];
    }
    return sum;
}
bool detectCycle(int* pre,int n, vector<int> *cycle){
    bool on_stack[n];
    bool visited[n];
    std::fill(visited, visited + n, false);
    std::fill(on_stack, on_stack + n, false);
    for (int i = 0; i < n; ++i) {
        if (!visited[i]) {
            for (int j = i; j != -1; j = pre[j]) {
                if (!visited[j]) {
                    visited[j] = true;
                    (*cycle).push_back(j);
                    on_stack[j] = true;
                } else {
                    if (on_stack[j]) {
                        //look if cycle remains correct
                        auto it = find((*cycle).begin(),(*cycle).end(),j);
                        (*cycle).erase((*cycle).begin(),it);
                        (*cycle).push_back(j);// to form the cycle
                        return true;
                    }
                    break;
                }
            }
            for (int j : *cycle)
                on_stack[j] = false;
            (*cycle).clear();
        }
    }
    return false;
}



int negativeCycleDetection(vector<unordered_map<int,int>> * residualGraph){
    int n = (*residualGraph).size();
    int dis[n];
    int len[n];
    int pre[n];
    bool in_queue[n];
    std::fill(dis, dis + n, 0);
    std::fill(len, len + n, 0);
    std::fill(pre, pre + n, -1);
    std::fill(in_queue, in_queue + n, true);
    queue<int> que;
    for (int i = 0; i < n; ++i)
        que.push(i);
    int iter = 0;
    while(!que.empty()){
        int u = que.front();
        que.pop();
        in_queue[u] = false;
        //string vName = "v" + to_string(uID);
        for(auto to : (*residualGraph)[u]){
            int v = to.first;
            if(dis[u]+to.second < dis[v]){
                pre[v] = u;
                len[v] = len[u] + 1;
                dis[v] = dis[u] + to.second;
                if(len[v] == n || ++iter == n){
                    iter = 0;
                    vector<int>cycle;
                    stack<string>sCycle;
                    if(detectCycle(pre,n,&cycle)){
                        //detect sycle should return true if there is cycle and should fill cycle stack
                        int cost = sum(&cycle,*residualGraph);
                        // for(auto v : cycle){
                        //     sCycle.push(v);
                        // }
                        // AugmentingPath *au = new AugmentingPath(1,0,sCycle);
                        // augment(residualGraph,*au);
                        //printGraph(residualGraph);
//                        for(auto v : cycle){
//                            cout << v << endl;
//                        }
                        return cost;
                    }
                }
                if(!in_queue[v]){
                    que.push(v);
                    in_queue[v] = true;
                }
            }
        }
    }
    return 0;
}



// Finds an augmenting path from source to target in the residual graph. Uses DFS to find the path.
// AugmentingPath findAugmentingPath(unordered_map<string, Vertex>& residualGraph, string sourceName, string targetName, string notReached="X") {
//     // Will be used to construct the path from source to target, when the target is reached.
//     unordered_map<string, string> reachedFrom;//predecessors
//     for (auto v: residualGraph) {
//         reachedFrom[v.first] = notReached;
//     }
//     int pre[];
//     std::fill(pre, pre + n, -1);

//     // This block is the DFS implementation and track where each vertex is reached from.
//     stack<string> dfsStack;
//     dfsStack.push(sourceName);
//     while(!dfsStack.empty() && reachedFrom[targetName] == notReached) {
//         string currentVertex = dfsStack.top();
//         dfsStack.pop();
//         for (auto edge: residualGraph[currentVertex].adjacencyMap) {
//             string toName = edge.first;
//             if (reachedFrom[toName] == notReached) {
//                 reachedFrom[toName] = currentVertex;
//                 dfsStack.push(toName);
//             }
//         }
//     }

//     // If the target is not reached, return an empty path.
//     if (reachedFrom[targetName] == notReached) {
//         return AugmentingPath();
//     }

//     // Backtrack from the target to construct the path. Use a stack to store the path for "LIFO" property.
//     // Also find the minimum capacity (bottleneck) along the path to decide the amount of flow.
//     int minCapacity = numeric_limits<int>::max();
//     int pathCost = 0;
//     stack<string> path;
//     string childName = targetName;
//     while (childName != sourceName) {
//         path.push(childName);
//         string parentName = reachedFrom[childName];
//         minCapacity = min(minCapacity, residualGraph[parentName].adjacencyMap[childName].first);
//         pathCost += residualGraph[parentName].adjacencyMap[childName].second;
//         childName = parentName;
//     }
//     path.push(sourceName);
//     return AugmentingPath(minCapacity, pathCost, path);
// }

// Update the residual graph to realize the augmentation. The path is stored as a path, so the first popped element is the source.
// void augment(unordered_map<string, Vertex>& residualGraph, AugmentingPath augmentingPath) {
//     int amount = augmentingPath.amount;//flow
//     //int cost = augmentingPath.totalCost;
//     stack<string> path = augmentingPath.path;
//     //cout << "Augmenting " << amount << " along ";
//     while (path.size() > 1) {
//         string fromVertex = path.top();
//         path.pop();
//         string toVertex = path.top();

//         // Decrease the forward capacity and delete the edge if the remaining capacity is 0.
//         int fwdCapacity = residualGraph[fromVertex].adjacencyMap[toVertex].first;
//         int fwdCost = residualGraph[fromVertex].adjacencyMap[toVertex].second;
//         int remainingCapacity = fwdCapacity - amount;

//         //cout << fromVertex << " - ";

//         if (remainingCapacity == 0) {
//             residualGraph[fromVertex].adjacencyMap.erase(toVertex);
//         } else {
//             residualGraph[fromVertex].adjacencyMap[toVertex].first = remainingCapacity;
//         }

//         // Add a backward edge or increment its capacity if it already exists.
//         // Note: at(val) returns an exception if val is not in the unordered_map.
//         try {
//             // There is already an edge. Increment the capacity
//             residualGraph[toVertex].adjacencyMap.at(fromVertex).first += amount;
//             //should i manage cost in already existing edge situation
//         } catch (const std::out_of_range& oor) {
//             // Add an edge.
//             residualGraph[toVertex].adjacencyMap[fromVertex].first = amount;
//             residualGraph[toVertex].adjacencyMap.at(fromVertex).second = fwdCost * (-1) ;
//         }
//     }
//     //cout << path.top() << endl;
// }


// Run Ford-Fulkerson algorithm to find the maximum flow in a graph.
// int fordFulkerson(unordered_map<string, Vertex>* graph, string sourceName, string targetName) {
//     string notReachedFlag = "X";
//     unordered_map<string, Vertex> residualGraph = *graph;
//     //cout << "Initial residual graph:\n" ;
//     //printGraph(residualGraph);

//     int maxFlow = 0, i = 0, totalCost=0;
//     // Augment until no path is left.
//     bool isAugmentingPathLeft = true;
//     while (isAugmentingPathLeft) {
//         AugmentingPath augmentingPath = findAugmentingPath(residualGraph, sourceName, targetName, notReachedFlag);

//         if (augmentingPath.amount > 0) {
//             augment(residualGraph, augmentingPath);
//             maxFlow = maxFlow + augmentingPath.amount;
//             totalCost += augmentingPath.totalCost;
//             i++;
//             //cout << "cost" << totalCost << endl;
//             // << "Residual graph at iteration: " << i << endl;
//             //printGraph(residualGraph);
//         } else {
//             isAugmentingPathLeft = false;
//         }
//     }
//     int neg = 0 ;
//     int add=2;

//     while(add!=0){
//         add = negativeCycleDetection(residualGraph);
//         neg+=add;
//     }
//     totalCost += neg;
//     cout << totalCost << endl;
//     return totalCost;
// }
int main(int argc, char* argv[]) {

    string infile_name = argv[1];//READ INPUT FILE NAME
    string outfile_name = argv[2];//READ OUTPUT FILE NAME

    //OPEN INPUT AND OUTPUT FILE STREAMS
    ifstream infile;
    infile.open(infile_name);
    ofstream outfile;
    outfile.open(outfile_name);

    int numOfCases;
    infile >> numOfCases;
    vector<vector<unordered_map<int,int>>> graphs;
    //vector<unordered_map<string, Vertex>> graphs;
    int ind = 0;
    vector<int> graphcosts(numOfCases);
    while(numOfCases > 0){
        int numOfCables;
        infile >> numOfCables;
        //graphSizes.push(numOfCables);
        //unordered_map<string, Vertex> graph;
        vector<unordered_map<int,int>> graph(2*numOfCables);
        // for (int i = 0; i < (2*numOfCables); ++i)
        // {
        //     string vertexName = "v" + to_string(i);
        //     graph[vertexName] = Vertex(vertexName,i);
        // }
        for (int r = 0; r < numOfCables; r++){
            for (int c = numOfCables; c < 2 * numOfCables; c++){
                int voltage;
                infile >> voltage;
                //cout << c << endl;
                //string from = "v" + to_string(r);
                // graph[from] = Vertex(from);
                //string to = "v"+to_string(c);
                //graph[from].addEdge(to,make_pair(1,voltage*(-1)));
                if(c == r+numOfCables){
                    graph[c][r] = voltage;//directly connecting edges parallel
                    graphcosts[ind]+=voltage*(-1);
                    //totalCost+=voltage*(-1);
                }
                else{
                    graph[r][c] = voltage*(-1);
                }
                //cout << from << " " << to << endl;
            }
        }
        // string s = "v" + to_string(2*numOfCables);
        // string t = "v" + to_string((2*numOfCables)+1);
        // graph[s] = Vertex(s,2*numOfCables);
        // graph[t] = Vertex(t,(2*numOfCables)+1);
        // int s = 2*numOfCables;
        // int t = (2*numOfCables)+1;
        // for (int i = 0; i < numOfCables; ++i)
        // {
        //     //graph[s].addEdge("v"+to_string(i),make_pair(1,0));
        // 	graph[s][i] = 0;
        //     //graph["v"+to_string(numOfCables+i)].addEdge(t,make_pair(1,0));
        // 	graph[numOfCables+i][t] = 0;
        // }
        graphs.push_back(graph);
        ind++;
        numOfCases--;
    }
    int i =0;
    for(auto g : graphs){
        // int gsize = graphSizes.front();
        // string s = "v" + to_string((2*gsize));
        // string t = "v" + to_string((2*gsize)+1);
        // graphSizes.pop();
        //int minCost = (-1)*fordFulkerson(&g, s, t);
        int neg = 0 ;
        int add=2;
        while(add!=0){
            add = negativeCycleDetection(&g);
            neg+=add;
        }
        graphcosts[i]+=neg;
        //totalCost += neg;
        //minCost += negativeCycleDetection(g);
        outfile << graphcosts[i] * -1 << endl;
        //cout << "Maximum flow is: " << fordFulkerson(&g, "s", "t") << endl;
        i++;
    }
    //  for(auto g : graphs){
    //  	for(auto from : g){
    // for(auto to : from){
    // 	cout << to.first << to.second <<  endl;
    // }
    //  	}
    //  }
    //cout << "Maximum flow is: " << fordFulkerson(graphs[0], "s", "t") << endl;
    // string v = "v" + to_string(0);
    // cout << v << endl;
    // for(auto g : graphs){
    // 	// vector<Edge*> g = graphs.front();
    // 	// graphs.pop();
    // 	for(auto edge : g){
    // 		cout << edge->from << edge->to << edge->cost << endl;
    // 	}
    // }
    // for(auto g : graphs){
    // 	// vector<Edge*> g = graphs.front();
    // 	// graphs.pop();
    // 	int gSize = graphSizes.front();
    // 	graphSizes.pop();
    // 	cout << gSize << endl;
    // 	int * distance = new int[gSize];
    // 	int * predecessor = new int[gSize];
    // 	for (int v = 0; v < gSize; ++v)
    // 	{
    // 		distance[v] = 10000;
    // 		predecessor[v] = -1;
    // 	}
    // 	distance[gSize-2] = 0;
    // 	for (int i = 0; i < gSize-1; ++i){
    // 		for (int i = 0; i < g.size(); ++i){
    // 			Edge* edge = g[i];
    // 			if(edge->capacity!=0){
    // 				if(distance[edge->from] + 1 < distance[edge->to]){
    // 					cout  << " edge" << i << edge->capacity << endl;
    // 					distance[edge->to] = distance[edge->from] + 1;
    // 					predecessor[edge->to] = edge->from;
    // 					edge->capacity = 0;
    // 				}
    // 			}
    // 		}
    // 	}
    // 	for(auto edge : g){
    // 		cout << edge->capacity <<endl;
    // 	}
    // 	for (int i = 0; i < gSize; ++i){
    // 		cout << "predecessor of " << i << " = " << predecessor[i] << endl;
    // 	}
    // }

    return 0;
}
