#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <queue>
#include <algorithm>
#include <unordered_map>
#include <string>
#include <stack>
using namespace std;

int sum(vector<int> *cycle,vector<unordered_map<int,int>> & residualGraph){
    int firstV = -1;
    int secondV = -1;
    int sum = 0;

    for (int i = (*cycle).size()-1; i >= 0; i--) {
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

// Modified the code from https://konaeakira.github.io/posts/using-the-shortest-path-faster-algorithm-to-find-negative-cycles.html
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


//  Modified the code from https://konaeakira.github.io/posts/using-the-shortest-path-faster-algorithm-to-find-negative-cycles.html
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
    int ind = 0;
    vector<int> graphcosts(numOfCases);
    while(numOfCases > 0){
        int numOfCables;
        infile >> numOfCables;
        vector<unordered_map<int,int>> graph(2*numOfCables);
        for (int r = 0; r < numOfCables; r++){
            for (int c = numOfCables; c < 2 * numOfCables; c++){
                int voltage;
                infile >> voltage;
                if(c == r+numOfCables){
                    graph[c][r] = voltage;//directly connecting edges parallel
                    graphcosts[ind]+=voltage*(-1);
                }
                else{
                    graph[r][c] = voltage*(-1);
                }
            }
        }
      
        graphs.push_back(graph);
        ind++;
        numOfCases--;
    }
    int i =0;
    for(auto g : graphs){
        int neg = 0 ;
        int add=2;
        while(add!=0){
            add = negativeCycleDetection(&g);
            neg+=add;
        }
        graphcosts[i]+=neg;
        outfile << graphcosts[i] * -1 << endl;
        i++;
    }
    

    return 0;
}
