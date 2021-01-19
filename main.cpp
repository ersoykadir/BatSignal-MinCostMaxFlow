#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <queue>
#include <unordered_map>
#include <string>
#include <stack>
#include "utils.cpp"
using namespace std;

bool stackContains(stack<int> st,int vID){
	if(st.size()==0) return false;
	while(!st.empty()){
		int v = st.top();
		st.pop();
		if(v == vID) return true;
	}
	return false;
}

int sum(vector <int> * cycle,unordered_map<string, Vertex>& residualGraph){
	auto iter = (*cycle).begin();
	int firstV,secondV,sum;
	while(iter!=(*cycle).end()){
		firstV = (*iter);
		secondV = *iter++;
		for(auto v : residualGraph){
			if(v.second.ID == firstV){
				for(auto to : v.second.adjacencyMap){
					if(secondV == residualGraph[to.first].ID){
						sum+=to.second.second;
					}
				}
			}
		}
	}
	return sum;
}

void negCycleCalc(int* pre, int vID,vector<int> * cycle){
	//int v = vID;
	stack<int> negCycleStack;
	while(!stackContains(negCycleStack,vID)){
		negCycleStack.push(vID);
		vID = pre[vID];
	}
	
	(*cycle).push_back(vID);
	while(negCycleStack.top() != vID){
		(*cycle).push_back(negCycleStack.top());
		negCycleStack.pop();
	}
	(*cycle).push_back(vID);
}

bool queContains(queue<int> que,int vID){// without using &que copyinh original que and working on copy i hope
	//queue<int> mockQ;
	if(que.size()==0) return false;
	while(!que.empty()){
		int v = que.front();
		que.pop();
		if(v == vID) return true;
	}
	return false;
}

int negativeCycleDetection(unordered_map<string, Vertex>& residualGraph){
	int * dis = new int [residualGraph.size()];
	int * len = new int [residualGraph.size()];
	int * pre = new int [residualGraph.size()];
	queue<int> que;
	for (auto v : residualGraph)
	{
		dis[v.second.ID] = 0;
		len[v.second.ID] = 0;
		pre[v.second.ID] = -1;
		que.push(v.second.ID);
	}
	while(!que.empty()){
		int uID = que.front();
		que.pop();
		string vName = "v" + to_string(uID);
		for(auto to : residualGraph[vName].adjacencyMap){
			int vID = residualGraph[to.first].ID;
			if(dis[uID]+to.second.second < dis[vID]){
				pre[vID] = uID;
				len[vID] = len[uID] + 1;
                if(len[vID] == residualGraph.size()){
                	vector<int> * cycle;
                	negCycleCalc(pre,vID,cycle);
                	return sum(cycle,residualGraph);
                }
                dis[vID] = dis[uID] + to.second.second;
                if(!queContains(que,vID)){
                	que.push(vID);
                }
			}
		}
	}
	return 0;
}



// Finds an augmenting path from source to target in the residual graph. Uses DFS to find the path.
AugmentingPath findAugmentingPath(unordered_map<string, Vertex>& residualGraph, string sourceName, string targetName, string notReached="X") {
	// Will be used to construct the path from source to target, when the target is reached.
	unordered_map<string, string> reachedFrom;//predecessors
	for (auto v: residualGraph) {
		reachedFrom[v.first] = notReached;
	}

	// This block is the DFS implementation and track where each vertex is reached from.
    stack<string> dfsStack;
    dfsStack.push(sourceName);
    while(!dfsStack.empty() && reachedFrom[targetName] == notReached) {
        string currentVertex = dfsStack.top();
        dfsStack.pop();
        for (auto edge: residualGraph[currentVertex].adjacencyMap) {
        	string toName = edge.first;
            if (reachedFrom[toName] == notReached) {
                reachedFrom[toName] = currentVertex;
                dfsStack.push(toName);
            }
        }
    }

    // If the target is not reached, return an empty path.
    if (reachedFrom[targetName] == notReached) {
    	return AugmentingPath();
    }

    // Backtrack from the target to construct the path. Use a stack to store the path for "LIFO" property.
    // Also find the minimum capacity (bottleneck) along the path to decide the amount of flow.
    int minCapacity = numeric_limits<int>::max();
    int pathCost = 0;
	stack<string> path;
	string childName = targetName;
	while (childName != sourceName) {
		path.push(childName);
		string parentName = reachedFrom[childName];
		minCapacity = min(minCapacity, residualGraph[parentName].adjacencyMap[childName].first);
		pathCost += residualGraph[parentName].adjacencyMap[childName].second;
		childName = parentName;
	}
	path.push(sourceName);
	return AugmentingPath(minCapacity, pathCost, path);
}

// Update the residual graph to realize the augmentation. The path is stored as a path, so the first popped element is the source.
void augment(unordered_map<string, Vertex>& residualGraph, AugmentingPath augmentingPath) {
	int amount = augmentingPath.amount;//flow
	//int cost = augmentingPath.totalCost;
	stack<string> path = augmentingPath.path;
	cout << "Augmenting " << amount << " along ";
	while (path.size() > 1) {
		string fromVertex = path.top();
		path.pop();
		string toVertex = path.top();

		// Decrease the forward capacity and delete the edge if the remaining capacity is 0.
		int fwdCapacity = residualGraph[fromVertex].adjacencyMap[toVertex].first;
		int fwdCost = residualGraph[fromVertex].adjacencyMap[toVertex].second;
		int remainingCapacity = fwdCapacity - amount;

		cout << fromVertex << " - ";

		if (remainingCapacity == 0) {
			residualGraph[fromVertex].adjacencyMap.erase(toVertex);
		} else {
			residualGraph[fromVertex].adjacencyMap[toVertex].first = remainingCapacity;
		}

		// Add a backward edge or increment its capacity if it already exists. 
		// Note: at(val) returns an exception if val is not in the unordered_map.
		try {
			// There is already an edge. Increment the capacity
			residualGraph[toVertex].adjacencyMap.at(fromVertex).first += amount;
			//should i manage cost in already existing edge situation
		} catch (const std::out_of_range& oor) {
			// Add an edge.
			residualGraph[toVertex].adjacencyMap[fromVertex].first = amount;
			residualGraph[toVertex].adjacencyMap.at(fromVertex).second = fwdCost * (-1);
		}
	}
	cout << path.top() << endl;
}


// Run Ford-Fulkerson algorithm to find the maximum flow in a graph.
int fordFulkerson(unordered_map<string, Vertex>* graph, string sourceName, string targetName) {
	string notReachedFlag = "X";
	unordered_map<string, Vertex> residualGraph = *graph;
	cout << "Initial residual graph:\n" ;
	printGraph(residualGraph);

	int maxFlow = 0, i = 0, totalCost=0;
	// Augment until no path is left.
	bool isAugmentingPathLeft = true;
	while (isAugmentingPathLeft) {
		AugmentingPath augmentingPath = findAugmentingPath(residualGraph, sourceName, targetName, notReachedFlag);

		if (augmentingPath.amount > 0) {
			augment(residualGraph, augmentingPath);
			maxFlow = maxFlow + augmentingPath.amount;
			totalCost += augmentingPath.totalCost;
			i++;
			cout << "cost" << totalCost << endl;
			cout << "Residual graph at iteration: " << i << endl;
			printGraph(residualGraph);
		} else {
			isAugmentingPathLeft = false;			
		}
	} 
	totalCost += negativeCycleDetection(residualGraph);
		//cout << minCost << endl;
	return totalCost;
}
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
	vector<unordered_map<string, Vertex>> graphs;
	queue<int> graphSizes;
	while(numOfCases > 0){
		int numOfCables;
		infile >> numOfCables;
		graphSizes.push((2*numOfCables)+2);
		unordered_map<string, Vertex> graph;
		for (int i = 0; i < (2*numOfCables); ++i)
		{
			string vertexName = "v" + to_string(i);
			graph[vertexName] = Vertex(vertexName,i);
		}
		for (int r = 0; r < numOfCables; r++){
			for (int c = numOfCables; c < 2 * numOfCables; c++){
				int voltage;
				infile >> voltage;
				//cout << c << endl;
				string from = "v" + to_string(r);
				// graph[from] = Vertex(from);
				string to = "v"+to_string(c);
				graph[from].addEdge(to,make_pair(1,voltage*(-1)));
				cout << from << " " << to << endl;
			}
		}
		graph["s"] = Vertex("s",2*numOfCables);
		graph["t"] = Vertex("t",2*numOfCables+1);
		for (int i = 0; i < numOfCables; ++i)
		{
			graph["s"].addEdge("v"+to_string(i),make_pair(1,0));
			graph["v"+to_string(numOfCables+i)].addEdge("t",make_pair(1,0));
		}
		graphs.push_back(graph);
		numOfCases--;
	}
	for(auto g : graphs){
		int minCost = fordFulkerson(&g, "s", "t");
		//minCost += negativeCycleDetection(g);
		//cout << minCost << endl;
		//cout << "Maximum flow is: " << fordFulkerson(&g, "s", "t") << endl;
		
	}
	// for(auto g : graphs){
	// 	for(auto v : g){
	// 		if(v.first == "v"+to_string(0)){
	// 			for(auto t : v.second.adjacencyMap){
	// 				cout << t.first << endl;
	// 			}
	// 		}
	// 	}
	// }
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
