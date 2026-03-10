#include<iostream>
#include<fstream>
#include<stdio.h>
#include<assert.h>
#include<vector>
#include<time.h>
#include <algorithm>
#include "queue.h"
//#include "nonUniformArray.h"
//#include "nonUniformArraySmaller.h"
using namespace std;

static clock_t timeAlg = 0;

static double zeta; // the large shortest distance of two unreachable nodes

// for DFS
static bool visited[MAX_NODES]; // whether it has been visited before

// for shortest distance computation
static  int distances[MAX_NODES]; // distances to a source node
static int numShortestPathsBefore[MAX_NODES];
static int numShortestPathsAfter[MAX_NODES];

static bool isFoundforGBC[MAX_NODES];
static double betweenCentralities[MAX_NODES];
static double currrentBetweenCentralities[MAX_NODES];

const size_t numSamples = 500000;
//const size_t maxNodesEachSample = MAX_ROWS;
//static int sampleShortestPathsArray[numSamples][maxNodesEachSample];

static int nodesInEachSampleSetS[numSamples];
int* sampleShortestPathsArraySetS[numSamples];
int totalSampledShortestPathsSetS = 0;

static int nodesInEachSampleSetT[numSamples];
int* sampleShortestPathsArraySetT[numSamples];
int totalSampledShortestPathsSetT = 0;

//static Queue  bfsQueue; // 

// for dominator algorithm: see paper `A fast algorithm for finding domnators in a flow graph'
static int reachableSet[MAX_NODES]; // the set of reachable nodes from some source node
static int numReachNodes; // number of Reachable nodes

static nonUniformArray<int> propagationGraph;
static nonUniformArray<int> predecessorList;


const size_t maxSampledSources = 12800;
nonUniformArraySmaller<int>* sourcesSuccessorsAdjList[maxSampledSources];
nonUniformArraySmaller<int>* shortestPathsAdjLists;
int* sourcesNumShortestPathsBefore[maxSampledSources];
int* oneSourceNumShortestPathsBefore;
int* sourcesNumShortestPathsAfter[maxSampledSources];
int* oneSourceNumShortestPathsAfter;
int* sourcesIn2OutIndex[maxSampledSources];

nonUniformArraySmaller<int>* sourcesPredecessorsAdjList[maxSampledSources];
nonUniformArraySmaller<int>* shortestPathsPredecessorAdjLists;
double* sourcesGBC[maxSampledSources];
double* oneSourceGBC;
int* oneSourceDecreasedNumShortestPaths;

static bool visitedBackwards[MAX_NODES];
static  int distancesBackwards[MAX_NODES];
static int numShortestPathsBackward[MAX_NODES]; 

struct PAIR{
	int s;
	int t;
	PAIR(int src, int sink) { s = src; t = sink; }
	PAIR() { s = t = 0; }
};

struct EdgeWithWeight {
	int u;
	int v;
	double weight;
	int direction; // 1: forward, 2: backward
	EdgeWithWeight(int startNode, int endNode, int di, double edgeWeight)
	{
		u = startNode;
		v = endNode; 
		direction = di;
		weight = edgeWeight;
	}
};

// another, fast way to represent the buckets
struct LinkNode{ 
	int node;
	int nextNodeIndex; // the index of next node in the list. -1 indicates that it is the last element
};

static int visitNumber = 0;

static bool isReachable[MAX_NODES];

static nonUniformArray<int> parentSCCNeighbors;
static nonUniformArray<int> sccNeighborGroups;
static nonUniformArray<int> sccNeighborSizes;


static int reachSCCsSet[MAX_NODES];
static int numReachSCCs;

static int globalDegrees[MAX_NODES];

static double sumEachSCC[MAX_NODES];

struct Node;
static vector<Node> sortSNs;

static nonUniformArray<int> oneSCC;
static nonUniformArray<int> SCCs;
static nonUniformArray<int> sccGraph;
static nonUniformArray<int> sccTransposeGraph;
static nonUniformArray<int> sccUnderlyingGraph;


// for relabelling the propagation graph 
static int outIndex2InIndex[MAX_NODES];
static int inIndex2OutIndex[MAX_NODES];

class GroupBetweenesssCentrality;

struct Node
{
	int id;
	double key;
	Node(int a, double b){ id = a; key = b; }
	Node(){}
	bool operator>(Node &e)
	{
		if (key > e.key) return true;
		else					return false;
	}
	bool operator>=(Node &e)
	{
		if (key > e.key) return true;
		else					return false;
	}
	bool operator<(Node &e)
	{
		if (key < e.key) return true;
		else					return false;
	}
	bool operator<=(Node &e)
	{
		if (key < e.key) return true;
		else					return false;
	}
};

inline bool compLess(Node& a, Node& b) // operator <
{
	return a.key < b.key;
}

inline bool comp(Node & a, Node & b) // operator >
{
	return a.key > b.key;
}

//static int Mode = 9973;
//static int Mode = 5499979;
inline double Uniform()
{
	//return (rand() % Mode) / (double)Mode;
	return double(rand()) / RAND_MAX;
}

int Uniform(int a, int b)
{
	return  floor(double(rand()) / RAND_MAX * (b - a) + a);
}

inline void swap(int &u, int &v)
{
	int tmp = u;
	u = v;
	v = tmp;
}

struct Edge{
	int u, v;
	Edge(int a, int b){ u = a; v = b; }
	Edge(){}
};

class Graph{
public:
	friend class GroupBetweenesssCentrality;
	friend class HAM;
public:
	void readGraph(const char *graphFile, bool isDirected);
	void generateEdgeNodeWeights(double minEdgeWeight, double maxEdgeWeight, double minNodeWeight, double maxNodeWeight);
	void generateEdgeNodeWeightsCaseStudy(double minNodeWeight, double maxNodeWeight);
	void printGraph();
	void printGraphTopology();
	int findAllCCsDFS(nonUniformArray<int> &allCCs); // find all connected components
	int findAllCCsBFS(nonUniformArray<int> &allCCs); // find all connected components
	void findAPsAndLowerBounds();
	void cleanGraph(const char *graphFile);
	Graph();
	Graph(const Graph &g);
	Graph & operator=(const Graph &g);
	~Graph();

public:
// for shortest distance
public:
	int numOfEdges();
	int numOfNodes(){ return n; }

protected:
	
	void rankByBetweennessCentrality(vector<int>& rankedNodes);
	double rankByGroupBCWithSampleSources(vector<int> &sourceNodes, vector<int>& foundNodes, int K);
	//double rankByGroupBCWithSampleSourcesFast(vector<int>& sourceNodes, vector<int>& foundNodes, int K);
	
	double rankByGroupBCWithSampleSourcesMemoryFriendly(vector<int>& sourceNodes, vector<int>& foundNodes, int K, vector<double>& foundTimes);

	double rankByGroupBCWithSampleSourcesMemoryFriendlyBFS(vector<int>& sourceNodes, vector<int>& foundNodes, int K, vector<double>& foundTimes);

	//double rankByGroupBCWithSampleSourcesMemoryFriendlyDFSPostOrder(vector<int>& sourceNodes, vector<int>& foundNodes, int K, vector<double>& foundTimes);
	//void findReachNodeByPostOrderDFS(int u, nonUniformArray<int>& adjList, vector<int>& nodesInOrder);

	double rankByGroupBCWithSampleSourcesFastPlus(vector<int>& sourceNodes, vector<int>& foundNodes, int K, vector<double> &  foundTimes);
	void recursiveCalGBC(int v);

	
	double rankByGroupBCWithSampleShortestPaths(vector<PAIR>& randomPairs, vector<int>& foundNodes, int K, 
		vector<double>& foundTimes);
	double rankByGroupBCWithSampleShortestPathsBiBFS(vector<PAIR>& randomPairs, vector<int>& foundNodes, int K,
		vector<double>& foundTimes);

	double rankByGroupBCWithSampleShortestPathsBiBFSunbiasedKDD23(vector<PAIR>& randomPairs, vector<int>& foundNodes, int K,
		vector<double>& foundTimes);
	double evaluateGBCbySampleShoretstPaths(vector<int>& foundNodes, vector<PAIR> &SamplePairs);
	
	double rankByGroupBCWithSamplePairs(vector<PAIR>& randomPairs, vector<int>& foundNodes, int K, vector<double>& foundTimes);

	
	double groupBCevaluation(vector<int>& foundNodes, double  nodesFraction);
	double groupBC_BatchEvaluation(vector<int>& foundSets, int numEachSet, double & minGBC, double & maxGBC, double
		scaleDownFactor, const char* outFileName, int numSamples, int testTimes, double nodesFraction, bool isWrite);
	void groupBC_IncrementEvaluation(vector<int>& foundNodes, int Delta, double nodesFraction, double scaleDownFactor,
		vector<double>& gbcValues);
	void generateRandomPermulation(int n, vector<int>& randomizedSeq);

	double findTopKGBCbyKDD16(int K, double epsilon, double errorProb, vector<int>& foundNodes, int& sampledPaths);

	double findTopKGBCbyKDD23(int K, double epsilon, double errorProb, vector<int>& foundNodes, int& sampledPaths);

	double findTopKGBCAdaptiveSampling(int K, double epsilon, double errorProb, vector<int>& foundNodes, int& sampledPaths);

	void generateRandomPairs(int n, int L, vector<PAIR>& randomizedPairs);
protected:
	// single-source shortest distance
	void shortestDistances(int src, bool markSPT);
protected:
	void modifiedDFS(int u, nonUniformArray<int> & allCCs, int numNodesValid);
	inline void DFS(int u, int CCID, int & sizeOfCC);
	void clearGraph();
protected:
	nonUniformArray<int> adjLists;
	nonUniformArray<double> adjWeights;
	vector<double> nodesWeights;
	int n, m; // number of nodes and edges in the graph
	// whether a node is valid in the graph, as it may be identified as a structural hole
	// and removed
	bool *isValid;

public: //for DFS
	void initilizeVariables();
	void clearVariables();
	void allocateVariablesForDFS();
private: // for DFS
	int *numChildren; // number of chidren in the DFS tree
	int *numNodesSubtree; // number of nodes in the subtree rooted at each node in the DFS tree
	int *parent; // the parent of each node in the DFS tree
	bool *isAP; // whether the node is an Articulation point
	int *discoveredTime; // the discovered time of each vertex
	int *lowestTime; // the smallest discovered time of any neighbor of u’s descendants (through a back-edge)
	double *lowerBounds; // the lower bound for each node
	int time; // global time counter

protected:
	vector<int> sizeofEachCC; // the number of nodes in each connected component
	int *ccBelongto; // the cc id that each node belongs to

private:
	bool debug; //
};



int Graph::numOfEdges()
{
	int edges = 0;
	for (int i = 0; i < adjLists.size(); ++i)
		edges += adjLists.size(i);
	return edges / 2;
}




void Graph::generateRandomPermulation(int n, vector<int>& randomizedSeq)
{ // generate a random permulation of n numbers 0, 1, ..., n-1
	assert(n > 0);
	randomizedSeq.clear();
	randomizedSeq.reserve(n);
	int i;
	for (i = 0; i < n; ++i)
		randomizedSeq.push_back(i);

	int tmp;
	int index;
	double p;
	for (i = n - 1; i > 0; --i)
	{
		p = Uniform();
		index = int(floor(p * i)); // a random number between 0 and i
		// swap randomizedSeq[index] and randomizedSeq[i]
		tmp = randomizedSeq[index];
		randomizedSeq[index] = randomizedSeq[i];
		randomizedSeq[i] = tmp;
		//printf("swap %d th node and %d th node\n", index, i);
	}
	/*for (i = 0; i < n; ++i)
	{
		printf("%d\t", randomizedSeq[i]);
		if ((i + 1) % 10 == 0) printf("\n");
	}*/
}

void Graph::groupBC_IncrementEvaluation(vector<int>& foundNodes, int Delta, double nodesFraction, double scaleDownFactor,
	vector<double>& gbcValues)
{
	assert(foundNodes.size() > 0 && Delta > 0);
	int totalNodes = foundNodes.size();
	int numSets = totalNodes / Delta;
	assert(numSets * Delta == totalNodes);

	int i, j, k;
	bool debug = true;
	nonUniformArray<int>  successorsAdjLists; // record successors in the shortest paths
	nonUniformArray<int>  successorsAdjListsCompact; // record successors in the shortest paths
	vector<int> degrees; // degree of each node 
	degrees.reserve(n);
	for (i = 0; i < n; ++i) degrees.push_back(adjLists.size(i));
	successorsAdjLists.allocateMemory(degrees);
	successorsAdjListsCompact.allocateMemory(degrees);
	int src;
	int u, v, w;
	// the number of nodes that their shortest distances to src have been found
	int sz;
	int* adj;
	int startIndex;

	double totalGBC = 0;
	double* gbcEachSet = new double[numSets];
	for (i = 0; i < numSets; ++i) gbcEachSet[i] = 0;

	int newDistance;
	int reachableNodes;

	for (i = 0; i < n; ++i) isFoundforGBC[i] = false;

	vector<int> randomSeq;
	srand(2333);
	generateRandomPermulation(n, randomSeq);

	vector<int> sourceNodes;
	sourceNodes.reserve(n);
	assert(nodesFraction > 0 && nodesFraction <= 1);
	for (i = 0; i < nodesFraction * n; ++i)
		sourceNodes.push_back(randomSeq[i]);
	//printf("%d source nodes for evaluation: ", int(sourceNodes.size()));
	//for (j = 0; j < sourceNodes.size(); ++j) printf("%d,\t", sourceNodes[j]);
	//printf("\n\n");

	//printf("total nodes; %d, choosen nodes: %d\n", n, sourceNodes.size());
	double maxDiameter = 0;
	int totalSources = sourceNodes.size();
	clock_t st = clock();
	for (i = 0; i < totalSources; ++i)
	{
		if (i % 500 == 0) printf("total: %d, evaluation for %d th source, used time: %.3lf\n", totalSources, i+1, (clock()-st)/1000.0);
		src = sourceNodes[i];

		for (j = 0; j < n; ++j) successorsAdjLists.clear(j);
		for (j = 0; j < n; ++j) numShortestPathsBefore[j] = 0;
		numShortestPathsBefore[src] = 1;
		// find the shortest paths from each node 
		for (j = 0; j < n; ++j) visited[j] = false;
		visited[src] = true;
		distances[src] = 0;
		queue_clear(); // clear the queue
		queue_push(src);
		reachableNodes = 1;
		while (queue_size() > 0) // bfs
		{
			u = queue_pop();

			sz = adjLists.size(u);
			adj = &adjLists.access(u, 0);
			newDistance = distances[u] + 1;
			for (j = 0; j < sz; ++j)
			{
				v = adj[j];
				if (false == visited[v])
				{
					visited[v] = true;
					distances[v] = newDistance;
					queue_push(v);
					reachableNodes++;
					if (newDistance > maxDiameter) maxDiameter = newDistance;
				}
				if (distances[v] == newDistance) {
					numShortestPathsBefore[v] += numShortestPathsBefore[u];
					successorsAdjLists.push_back(u, v);
				}
			}
		}
		degrees.clear();
		for (j = 0; j < n; ++j) degrees.push_back(successorsAdjLists.size(j));
		successorsAdjListsCompact.allocateMemory(degrees);
		for (j = 0; j < n; ++j)
		{
			sz = successorsAdjLists.size(j);
			adj = &successorsAdjLists.access(j, 0);
			for (k = 0; k < sz; ++k)
				successorsAdjListsCompact.push_back(j, adj[k]);
		}

		for (j = 1; j <= numSets; ++j)
		{
			for (k = 0; k < j * Delta; ++k) isFoundforGBC[foundNodes[k]] = true;
			if (true == isFoundforGBC[src])
			{
				gbcEachSet[j - 1] += reachableNodes - 1;
				for (k = 0; k < j * Delta; ++k) isFoundforGBC[foundNodes[k]] = false;
				continue;
			}
			
			for (k = 0; k < n; ++k) numShortestPathsAfter[k] = 0;
			numShortestPathsAfter[src] = 1;
			queue_clear(); // clear the queue
			queue_push(src);
			for (k = 0; k < n; ++k) visited[k] = false;
			visited[src] = true;
			for (k = 0; k < j * Delta; ++k) visited[ foundNodes[k] ] = true;

			while (queue_size() > 0) // bfs
			{
				u = queue_pop();
				sz = successorsAdjListsCompact.size(u);
				adj = &successorsAdjListsCompact.access(u, 0);

				for (k = 0; k < sz; ++k)
				{
					v = adj[k];
					if (false == visited[v])
					{
						visited[v] = true;
						queue_push(v);
						numShortestPathsAfter[v] += numShortestPathsAfter[u];
					}
					else if (false == isFoundforGBC[v])
						numShortestPathsAfter[v] += numShortestPathsAfter[u];
				}
			}

			for (k = 0; k < n; ++k)
				if (numShortestPathsBefore[k] > 0)
					gbcEachSet[j-1] += 1.0 - double(numShortestPathsAfter[k]) / numShortestPathsBefore[k];

			for (k = 0; k < j * Delta; ++k) isFoundforGBC[foundNodes[k]] = false;
		}// for (j = 1; j <= numSets; ++j)
	}//for (i = 0; i < totalSources; ++i)
	//printf("maxDiameter: %.0lf\n", maxDiameter);
	double scaleUpFactor = double(n) / sourceNodes.size();
	printf("scaleUpFactor: %.3lf\n", scaleUpFactor);
	for (i = 0; i < numSets; ++i)
		gbcEachSet[i] *= scaleUpFactor;
	double alpha = double(n) * (n - 1.0);
	for (i = 0; i < numSets; ++i)
		gbcEachSet[i] /= alpha;
	for (i = 0; i < numSets; ++i)
		gbcEachSet[i] *= scaleDownFactor;
	
	gbcValues.clear();
	gbcValues.reserve(numSets);
	for (i = 0; i < numSets; ++i)
		gbcValues.push_back(gbcEachSet[i]);

	for (i = 0; i < numSets; ++i) printf("top %d, gbc: %.6lf\n", (i+1)*Delta, gbcEachSet[i]);
	delete[]gbcEachSet;
}
double Graph::groupBC_BatchEvaluation(vector<int>& foundSets, int numEachSet, double& minGBC, double& maxGBC,
	double scaleDownFactor, const char *outFileName, int numSamples, int testTimes, double nodesFraction, bool isWrite)
{
	int K = numEachSet;
	int numSets = foundSets.size() / K;
	assert(numSets * K == foundSets.size());

	int i, j, k;
	bool debug = true;
	nonUniformArray<int>  successorsAdjLists; // record successors in the shortest paths
	nonUniformArray<int>  successorsAdjListsCompact; // record successors in the shortest paths
	vector<int> degrees; // degree of each node 
	degrees.reserve(n);
	for (i = 0; i < n; ++i) degrees.push_back(adjLists.size(i));
	successorsAdjLists.allocateMemory(degrees);
	successorsAdjListsCompact.allocateMemory(degrees);
	int src;
	int u, v, w;
	// the number of nodes that their shortest distances to src have been found
	int sz;
	int* adj;
	int startIndex;

	double totalGBC = 0;
	double* gbcEachSet = new double[numSets];
	for (i = 0; i < numSets; ++i) gbcEachSet[i] = 0;
	/*if (debug)
	{
		printf("num nodes: %d\n", n);
		startIndex = 0;
		for (i = 0; i < numSets; ++i)
		{
			printf("set %d: ", i + 1);
			for (j = 0; j < K; ++j, ++startIndex) printf("%d, ", foundSets[startIndex]);
			printf("\n");
		}
		printf("\n");
	}*/
	
	int newDistance;
	int reachableNodes;
	
	for (i = 0; i < n; ++i) isFoundforGBC[i] = false;

	vector<int> randomSeq;
	srand(9999+ nodesFraction * (n+1));
	generateRandomPermulation(n, randomSeq);

	vector<int> sourceNodes;
	sourceNodes.reserve(n);
	assert(nodesFraction > 0 && nodesFraction <= 1);
	for (i = 0; i < nodesFraction * n; ++i)
		sourceNodes.push_back(randomSeq[i]);

	printf("total nodes; %d, choosen nodes: %d\n", n, sourceNodes.size());
	double maxDiameter = 0;
	for (i = 0; i < sourceNodes.size(); ++i)
	{	if( i %1000 == 0) printf("evaluation for source %d\n", i);
		src = sourceNodes[i];
	
		for (j = 0; j < n; ++j) successorsAdjLists.clear(j);
		for (j = 0; j < n; ++j) numShortestPathsBefore[j] = 0;
		numShortestPathsBefore[src] = 1;
		// find the shortest paths from each node 
		for (j = 0; j < n; ++j) visited[j] = false;
		visited[src] = true;
		distances[src] = 0;
		queue_clear(); // clear the queue
		queue_push(src);
		reachableNodes = 1;
		while (queue_size() > 0) // bfs
		{
			u = queue_pop();

			sz = adjLists.size(u);
			adj = &adjLists.access(u, 0);
			newDistance = distances[u] + 1;
			for (j = 0; j < sz; ++j)
			{
				v = adj[j];
				if (false == visited[v])
				{
					visited[v] = true;
					distances[v] = newDistance;
					queue_push(v);
					reachableNodes++;
					if (newDistance > maxDiameter) maxDiameter = newDistance;
				}
				if (distances[v] == newDistance) {
					numShortestPathsBefore[v] += numShortestPathsBefore[u];
					successorsAdjLists.push_back(u, v);
				}
			}
		}
		degrees.clear();
		for (j = 0; j < n; ++j) degrees.push_back(successorsAdjLists.size(j));
		successorsAdjListsCompact.allocateMemory(degrees);
		for (j = 0; j < n; ++j)
		{
			sz = successorsAdjLists.size(j);
			adj = &successorsAdjLists.access(j, 0);
			for (k = 0; k < sz; ++k)
				successorsAdjListsCompact.push_back(j, adj[k]);
		}

		startIndex = 0;
		for (j = 0; j < numSets; ++j, startIndex += K)
		{
			for (k = 0; k < K; ++k) isFoundforGBC[foundSets[startIndex + k]] = true;
			if (true == isFoundforGBC[src])
			{
				gbcEachSet[j] += reachableNodes - 1;
				for (k = 0; k < K; ++k)
					isFoundforGBC[foundSets[startIndex + k]] = false;
				continue;
			}

			for (k = 0; k < n; ++k) numShortestPathsAfter[k] = 0;
			numShortestPathsAfter[src] = 1;
			queue_clear(); // clear the queue
			queue_push(src);
			for (k = 0; k < n; ++k) visited[k] = false;
			visited[src] = true;
			for (k = 0; k < K; ++k) visited[foundSets[startIndex + k]] = true;

			while (queue_size() > 0) // bfs
			{
				u = queue_pop();
				//sz = successorsAdjLists.size(u);
				//adj = &successorsAdjLists.access(u, 0);
				sz = successorsAdjListsCompact.size(u);
				adj = &successorsAdjListsCompact.access(u, 0);

				for (k = 0; k < sz; ++k)
				{
					v = adj[k];
					if (false == visited[v])
					{
						visited[v] = true;
						queue_push(v);
						numShortestPathsAfter[v] += numShortestPathsAfter[u];
					}
					else if (false == isFoundforGBC[v])
						numShortestPathsAfter[v] += numShortestPathsAfter[u];
				}
			}

			for (k = 0; k < n; ++k)
				if (numShortestPathsBefore[k] > 0)
					gbcEachSet[j] += 1.0 -  double(numShortestPathsAfter[k]) / numShortestPathsBefore[k];

			for (k = 0; k < K; ++k) 
				isFoundforGBC[foundSets[startIndex + k]] = false;
		}// for (j = 0; j < numSets; ++j, startIndex += K)
	}// end for(i = 0; i < n; ++i)
	//printf("maxDiameter: %.0lf\n", maxDiameter);
	double scaleUpFactor = double(n) / sourceNodes.size();
	//printf("scaleUpFactor: %.3lf\n", scaleUpFactor);
	for (i = 0; i < numSets; ++i)
	 gbcEachSet[i] *= scaleUpFactor;
	double alpha = n * (n - 1.0);
	for (i = 0; i < numSets; ++i)
		gbcEachSet[i] /= alpha;
	for (i = 0; i < numSets; ++i)
		gbcEachSet[i] *= scaleDownFactor;
	totalGBC = 0;
	for (i = 0; i < numSets; ++i) totalGBC += gbcEachSet[i];
	double avgGBC = totalGBC / numSets;

	minGBC = 1.0;
	for (i = 0; i < numSets; ++i)
		if (gbcEachSet[i] <  minGBC) minGBC = gbcEachSet[i];
	maxGBC = 0;
	for (i = 0; i < numSets; ++i)
		if (gbcEachSet[i] > maxGBC) maxGBC = gbcEachSet[i];

	if (true == isWrite) {
		const int scale = 120;
		double prob[scale];
		double Delta = 0.6 / scale;
		for (i = 0; i < scale; ++i) prob[i] = 0;
		for (i = 0; i < numSets; ++i)
			prob[(int)floor(gbcEachSet[i] / Delta)]++;
		for (i = 0; i < scale; ++i)
			prob[i] /= numSets;
		char FileName[100];
		sprintf_s(FileName, "%s-%d-K=%d-L=%d-probDistribution.txt", outFileName, testTimes, K, numSamples);
		ofstream out(FileName, ios::app);
		for (i = 0; i < scale; ++i)
			out << i * Delta << '\t' << 100 * prob[i] << '\n';
		out.close();

		sprintf_s(FileName, "%s-%d-K=%d-L=%d.txt", outFileName, testTimes, K, numSamples);
		out.open(FileName, ios::app);
		out << numSamples << '\t' << minGBC << '\t' << avgGBC << '\t' << maxGBC << '\n';
		out.close();
		//printf("probability distribution:\n");
	//for (i = 0; i < scale; ++i)printf("[%.3lf, %.3lf), prob: %.2lf%%\n", i* Delta, (i+1)* Delta, 100*prob[i]);
	/**/
	}

	//printf("Evaluation: total gbc: %.6lf\n", totalGBC);
	delete []gbcEachSet;
	return avgGBC;

}
double Graph::groupBCevaluation(vector<int>& foundNodes, double  nodesFraction)
{
	int K = foundNodes.size();

	int i, j, k, l;
	bool debug = false;
	nonUniformArray<int>  successorsAdjLists; // record successors in the shortest paths
	vector<int> degrees; // degree of each node 
	degrees.reserve(n);
	for (i = 0; i < n; ++i) degrees.push_back(adjLists.size(i));
	successorsAdjLists.allocateMemory(degrees);
	int src;
	int u, v, w;
	// the number of nodes that their shortest distances to src have been found
	int sz;
	int* adj;

	double totalGBC = 0;

	//printf("num nodes: %d\n", n);
	int newDistance;
	int reachableNodes;
	for (i = 0; i < n; ++i) isFoundforGBC[i] = false;
	for (i = 0; i < K; ++i) isFoundforGBC[ foundNodes[i] ] = true;
	if (debug)
	{
		printf("GBC, found nodes: ");
		for (i = 0; i < K; ++i) printf("%d,\t", foundNodes[i]);
		printf("\n");
	}

	vector<int> randomSeq;
	srand(0);
	generateRandomPermulation(n, randomSeq);

	vector<int> sourceNodes;
	sourceNodes.reserve(n);
	assert(nodesFraction > 0 && nodesFraction <= 1);
	for (i = 0; i < nodesFraction * n; ++i)
		sourceNodes.push_back(randomSeq[i]);

	for (i = 0; i < sourceNodes.size(); ++i)
	{
		//src = i;
		src = sourceNodes[i];
		for (j = 0; j < n; ++j) successorsAdjLists.clear(j);
		for (j = 0; j < n; ++j) numShortestPathsBefore[j] = 0;
		numShortestPathsBefore[src] = 1;
		// find the shortest paths from each node 
		for (j = 0; j < n; ++j) visited[j] = false;
		visited[src] = true;
		distances[src] = 0;
		queue_clear(); // clear the queue
		queue_push(src);
		reachableNodes = 1;
		while (queue_size() > 0) // bfs
		{
			u = queue_pop();

			sz = adjLists.size(u);
			adj = &adjLists.access(u, 0);
			newDistance = distances[u] + 1;
			for (j = 0; j < sz; ++j)
			{
				v = adj[j];
				if (false == visited[v])
				{
					visited[v] = true;
					distances[v] = newDistance;
					queue_push(v);
					reachableNodes++;
				}
				if (distances[v] == newDistance) {
					numShortestPathsBefore[v] += numShortestPathsBefore[u];
					successorsAdjLists.push_back(u, v);
				}
			}
		}
		if (debug && src == 0)
		{
			printf("from source: %d, before\n", src);
			successorsAdjLists.printArray();
			for (j = 0; j < n; ++j)
				printf("to node %d, numShortestBefore: %d\n", j, numShortestPathsBefore[j]);
		}
		if (true == isFoundforGBC[src])
		{
			totalGBC += reachableNodes - 1;
			continue;
		}
			
		for (j = 0; j < n; ++j) numShortestPathsAfter[j] = 0;
		numShortestPathsAfter[src] = 1;
		queue_clear(); // clear the queue
		queue_push(src);
		for (j = 0; j < n; ++j) visited[j] = false;
		visited[src] = true;
		for (j = 0; j < foundNodes.size(); ++j) visited[foundNodes[j]] = true;

		while (queue_size() > 0) // bfs
		{
			u = queue_pop();
			sz = successorsAdjLists.size(u);
			adj = &successorsAdjLists.access(u, 0);

			for (j = 0; j < sz; ++j)
			{
				v = adj[j];
				if (false == visited[v])
				{
					visited[v] = true;
					queue_push(v);
					numShortestPathsAfter[v] += numShortestPathsAfter[u];
				}
				else if (false == isFoundforGBC[v])
					numShortestPathsAfter[v] += numShortestPathsAfter[u];
			}
		}
		if (debug && src == 0)
		{
			printf("from source: %d, after\n", src);
			successorsAdjLists.printArray();
			for (j = 0; j < n; ++j)
				printf("to node %d, numShortestAfter: %d\n", j, numShortestPathsAfter[j]);
		}
		
		for (j = 0; j < n; ++j)
			if (numShortestPathsBefore[j] > 0)
				totalGBC += 1.0 - 1.0*numShortestPathsAfter[j] / numShortestPathsBefore[j];
	}// end for(i = 0; i < n; ++i)
		
	//printf("Evaluation: total gbc: %.6lf\n", totalGBC/(n * (n - 1)));
	double scaleUpFactor = double(n) / sourceNodes.size();
	scaleUpFactor *= scaleUpFactor;
	return totalGBC/(n*(n-1.0));
}


double Graph::rankByGroupBCWithSampleShortestPathsBiBFSunbiasedKDD23(vector<PAIR>& randomPairs, vector<int>& foundNodes, int K, vector<double>& foundTimes)
{ // find the sample with bidirectional BFS
	clock_t  st = clock();

	foundNodes.clear();
	foundNodes.reserve(K);
	foundTimes.clear();
	assert(K > 0 && K <= n);
	int L = randomPairs.size();
	assert(L <= numSamples);
	PAIR onePair;

	int i, j, k;
	int destID;
	bool debug = false;

	int sz;
	int* adj;
	int u, v;

	vector<int> degrees; // degree of each node 
	degrees.reserve(n);
	for (i = 0; i < n; ++i) degrees.push_back(adjLists.size(i));
 
	int** originalGraph = new int* [n];
	for (i = 0; i < n; ++i) originalGraph[i] = new int[degrees[i] + 1]; // some times degrees[i] is zero
	int* originalGraphSize = new int[n];
	for (i = 0; i < n; ++i) originalGraphSize[i] = 0;
	for (i = 0; i < n; ++i)
	{
		sz = adjLists.size(i);
		adj = &adjLists.access(i, 0);
		for (j = 0; j < sz; ++j) {
			v = adj[j];
			originalGraph[i][j] = v;
		}
		originalGraphSize[i] = sz;
	}

	int  **predecessorsAdjLists = new int *[n]; // record predecessors in the shortest paths
	for (i = 0; i < n; ++i)
		predecessorsAdjLists[i] = new int[degrees[i] + 1];
	int* predecessorsAdjListsSize = new int[n];

	int src, dst;
	// the number of nodes that their shortest distances to src have been found
	//int numDistancesFound = 0;
	

	double totalGBC = 0;
	double largestIncreasedBC;
	int largestNodeID;

	//printf("num nodes: %d\n", n);
	int newDistance;
	double p;
	double delta;
	double numEdegesExplord = 0;

	int* sampleShortestPathsArray[numSamples];
	const int maxNodesInPath = 100;
	for (i = 0; i < L; ++i) sampleShortestPathsArray[i] = new int[maxNodesInPath];

	degrees.clear();
	for (i = 0; i < n; ++i) degrees.push_back(0); // degrees in the transpose graph
	for (i = 0; i < n; ++i)
	{
		sz = adjLists.size(i);
		adj = &adjLists.access(i, 0);
		for (j = 0; j < sz; ++j)
			degrees[adj[j]]++;
	}
	int** transposeGraph = new int* [n];
	for (i = 0; i < n; ++i) transposeGraph[i] = new int[degrees[i] + 1];
	int* transposeGraphSize = new int[n];
	for (i = 0; i < n; ++i)transposeGraphSize[i] = 0;

	for (i = 0; i < n; ++i)
	{
		sz = adjLists.size(i);
		adj = &adjLists.access(i, 0);
		for (j = 0; j < sz; ++j) {
			v = adj[j];
			transposeGraph[ v ][ transposeGraphSize[v] ] = i;
			transposeGraphSize[v]++;
		}
	}
	int** predecessorAdjListsBackwards = new int* [n];
	for (i = 0; i < n; ++i) predecessorAdjListsBackwards[i] = new int[degrees[i]+1];
	int* predecessorAdjListsBackwardsSize = new int[n];

	//printf("original graph:\n"); adjLists.printArray();
	//printf("\n transpose graph:\n"); transposeGraph.printArray();

	
	clock_t ft = st;
	int middleNode;
	int tempNode;
	vector<int> randomSeq;
	randomSeq.reserve(100);

	NodeInTwoQuques nodeU, nodeV;
	int stopBFSdistanceForward, stopBFSdistanceBackward;
	vector<EdgeWithWeight> touchEdges;
	double totalWeight;
	double sumWeight;

	for (i = 0; i < L; ++i) // randomly generate L shortest paths
	{
		//if (i == 40) debug = true;
		//else debug = false;
		//debug = true;

		onePair = randomPairs[i];
		src = onePair.s;
		dst = onePair.t;

		/*if (i % 1000 == 0) {
			printf("L: %d, find the %d th sample, time used: %.2lf, delta time: %.2lf\n", L, i + 1, (clock() - st) / 1000.0,
				(clock() - ft) / 1000.0);
			ft = clock();
		}*/
		for (j = 0; j < n; ++j) predecessorsAdjListsSize[j] = 0; 
		for (j = 0; j < n; ++j) predecessorAdjListsBackwardsSize[j] = 0;

		// find the shortest paths from src 
		for (j = 0; j < n; ++j) visited[j] = false;
		visited[src] = true;
		for (j = 0; j < n; ++j) visitedBackwards[j] = false;
		visitedBackwards[dst] = true;
		for (j = 0; j < n; ++j) numShortestPathsBefore[j] = 0;
		numShortestPathsBefore[src] = 1;
		for (j = 0; j < n; ++j) numShortestPathsBackward[j] = 0;
		numShortestPathsBackward[dst] = 1;

		distances[src] = 0;
		distancesBackwards[dst] = 0;
		biQueue_clear(); // clear the queue
		biQueue_push(NodeInTwoQuques(src, 1));  //first queue
		biQueue_push(NodeInTwoQuques(dst, 2));  // second queue
		stopBFSdistanceForward = stopBFSdistanceBackward = n - 1;
		touchEdges.clear();
		if (debug) printf("\npair(%d, %d)******************\n", src, dst);

		while (biQueue_size() > 0) // bfs
		{
			nodeU = biQueue_pop();
			u = nodeU.v;
			if (debug) printf("exploring node %d ", u);

			if (1 == nodeU.qID) // forward search
			{
				if (debug) printf("forward \n");
				if (distances[u] >= stopBFSdistanceForward) continue;

				sz = originalGraphSize[u]; 
				adj = &originalGraph[u][0]; 
				newDistance = distances[u] + 1;

				generateRandomPermulation(sz, randomSeq);
				for (j = 0; j < sz; ++j)
				{
					v = adj[randomSeq[j]]; 
					v = adj[j];
					if (false == visited[v])
					{
						visited[v] = true;
						distances[v] = newDistance;
						biQueue_push(NodeInTwoQuques(v, 1));
						if (debug) printf("node %d is added to queue\n", v);
					}
					if (distances[v] == newDistance) {
						predecessorsAdjLists[v][predecessorsAdjListsSize[v]] = u;
						predecessorsAdjListsSize[v]++;
						numShortestPathsBefore[v] += numShortestPathsBefore[u];
					}

					if (true == visitedBackwards[v])
					{
						stopBFSdistanceForward = newDistance;
						stopBFSdistanceBackward = distancesBackwards[v];
						touchEdges.push_back(EdgeWithWeight(u, v, 1, 0));  // forward edge
						if (debug) printf("at node %d, found a middle node: %d\n", u, v);
					}
				}
				numEdegesExplord += sz;
			}
			else { // backward search
				if (debug) printf("backward \n");
				if (distancesBackwards[u] >= stopBFSdistanceBackward) continue;
				sz = transposeGraphSize[u]; 
				adj = &transposeGraph[u][0];
				newDistance = distancesBackwards[u] + 1;
				generateRandomPermulation(sz, randomSeq);
				for (j = 0; j < sz; ++j)
				{
					v = adj[randomSeq[j]];
					//v = adj[j];
					if (false == visitedBackwards[v])
					{
						visitedBackwards[v] = true;
						distancesBackwards[v] = newDistance;
						biQueue_push(NodeInTwoQuques(v, 2));
						if (debug) printf("node %d is added to queue\n", v);
					}
					if (distancesBackwards[v] == newDistance) {
						predecessorAdjListsBackwards[v][predecessorAdjListsBackwardsSize[v]] = u;
						predecessorAdjListsBackwardsSize[v]++;
						numShortestPathsBackward[v] += numShortestPathsBackward[u];
					}

					if (true == visited[v])
					{
						stopBFSdistanceForward = distances[v];
						stopBFSdistanceBackward = newDistance;
						touchEdges.push_back(EdgeWithWeight(u, v, 2, 0)); // backward edge
						if (debug) printf("at node %d, found a middle node: %d\n", u, v);
					}
				}
				numEdegesExplord += sz;
			}
		}
		destID = onePair.t;
		nodesInEachSampleSetS[i] = 0;
		if (touchEdges.size() == 0) {
			printf("pair: %d, %d, not reachable\n", src, destID);
			sampleShortestPathsArray[i][0] = src;
			nodesInEachSampleSetS[i] = 1;
			continue;
		}

		totalWeight = 0;
		for (j = 0; j < touchEdges.size(); ++j)
		{
			if (1 == touchEdges[j].direction) // forward edge
			{
				u = touchEdges[j].u;
				v = touchEdges[j].v;
			}
			else { // backward edge
				u = touchEdges[j].v;
				v = touchEdges[j].u;
			}
			touchEdges[j].weight = double(numShortestPathsBefore[u]) * numShortestPathsBackward[v];
			totalWeight += touchEdges[j].weight;
			if (debug) printf("touch edge: (%d, %d), weight: %.1lf\n", u, v, touchEdges[j].weight);
		}
		p = Uniform();
		if (debug) printf("p: %.3lf\n", p);
		p *= totalWeight;
		sumWeight = 0;

		for (j = 0; j < touchEdges.size(); ++j)
		{
			sumWeight += touchEdges[j].weight;
			if (sumWeight >= p) break;
		}
		assert(j != touchEdges.size());
		middleNode = touchEdges[j].v;
		if (debug) printf("chosen edge: (% d, % d), weight: % .1lf\n", touchEdges[j].u, touchEdges[j].v, touchEdges[j].weight);

		if (debug)
		{
			printf("pair: (%d, %d)\n predecessorsAdjLists:\n", src, onePair.t);
			//predecessorsAdjLists.printArray();
			printf("predecessorAdjListsBackwords:\n");
			//predecessorAdjListsBackwards.printArray();
			printf("middle node: %d\n", middleNode);
		}

		tempNode = middleNode;
		while (tempNode != src)
		{
			assert(tempNode >= 0 && tempNode < n);
			if (nodesInEachSampleSetS[i] >= maxNodesInPath) abort();
			sampleShortestPathsArraySetS[i][nodesInEachSampleSetS[i]] = tempNode;
			nodesInEachSampleSetS[i]++;
			p = Uniform();
			sz = predecessorsAdjListsSize[tempNode];
			if (sz == 0) {
				printf("error: pair: %d, %d\npredecessorsAdjLists:\n", onePair.s, onePair.t);
				//predecessorsAdjLists.printArray();
				printf("%d th sample\n", i + 1);
				abort();
			}
			delta = 1.0 / sz;
			j = int(floor(p / delta)); // round bin gaming
			if (j == sz) j--;
			assert(j >= 0 && j < sz);
			//if (debug) printf("---INSERT node %d\n", tempNode);

			tempNode = predecessorsAdjLists[tempNode][j];
		}
		if (nodesInEachSampleSetS[i] >= maxNodesInPath) abort();
		sampleShortestPathsArraySetS[i][nodesInEachSampleSetS[i]] = src;
		nodesInEachSampleSetS[i]++;
		//if (debug) printf("---INSERT node %d\n", src);

		tempNode = middleNode;
		if (tempNode == dst) continue;
		p = Uniform();
		sz = predecessorAdjListsBackwardsSize[tempNode];
		delta = 1.0 / sz;
		j = int(floor(p / delta)); // round bin gaming
		if (j == sz) j--;
		assert(j >= 0 && j < sz);
		tempNode = predecessorAdjListsBackwards[tempNode][j];

		while (tempNode != dst)
		{
			assert(tempNode >= 0 && tempNode < n);
			if (nodesInEachSampleSetS[i] >= maxNodesInPath) abort();
			sampleShortestPathsArraySetS[i][nodesInEachSampleSetS[i]] = tempNode;
			nodesInEachSampleSetS[i]++;
			p = Uniform();
			sz = predecessorAdjListsBackwardsSize[tempNode];
			if (sz == 0) {
				printf("error: pair: %d, %d \n predecessorAdjListsBackwards:\n", onePair.s, onePair.t);
				//predecessorAdjListsBackwards.printArray();
				printf("%d th sample\n", i + 1);
				abort();
			}
			delta = 1.0 / sz;
			j = int(floor(p / delta)); // round bin gaming
			if (j == sz) j--;
			//if (j < 0) j = 0;
			assert(j >= 0 && j < sz);
			//if (debug) printf("---INSERT node %d\n", tempNode);

			tempNode = predecessorAdjListsBackwards[tempNode][j];
		}
		if (nodesInEachSampleSetS[i] >= maxNodesInPath) abort();
		sampleShortestPathsArraySetS[i][nodesInEachSampleSetS[i]] = dst;
		nodesInEachSampleSetS[i]++;
		//if (debug) printf("---INSERT node %d\n", dst);
		if (debug)
		{
			printf("%d th sample: ", i + 1);
			for (j = 0; j < nodesInEachSampleSetS[i]; ++j)
				printf("%d,", sampleShortestPathsArraySetS[i][j]);
			printf(", middle Point: %d\n", middleNode);
		}
	}
	foundTimes.push_back((clock() - st) / 1000.0);
	numEdegesExplord /= L;
	//printf("bi BFS numEdegesExplord: %.0lf\n", numEdegesExplord);

	for (i = 0; i < n; ++i) isFoundforGBC[i] = false;
	for (k = 0; k < K; ++k)
	{
		//if (k % 10 == 0) printf("finding the %d node...\n", k + 1);
		for (i = 0; i < n; ++i) betweenCentralities[i] = 0;
		for (i = 0; i < L; ++i)
		{
			for (j = 0; j < nodesInEachSampleSetS[i]; ++j)
				if (isFoundforGBC[sampleShortestPathsArray[i][j]] == true) break;
			if (j != nodesInEachSampleSetS[i]) continue; // already covered
			for (j = 0; j < nodesInEachSampleSetS[i]; ++j)
				betweenCentralities[sampleShortestPathsArray[i][j]]++;
		}
		largestIncreasedBC = 0;
		for (i = 0; i < n; ++i)
		{
			if (betweenCentralities[i] > 0)
				assert(isFoundforGBC[i] == false);
			if (betweenCentralities[i] > largestIncreasedBC)
			{
				largestIncreasedBC = betweenCentralities[i];
				largestNodeID = i;
			}
		}
		if (largestIncreasedBC < 0.5) break;// no nodes are found
		isFoundforGBC[largestNodeID] = true;
		foundNodes.push_back(largestNodeID);
		totalGBC += largestIncreasedBC;

		if ((k + 1) % 10 == 0) foundTimes.push_back((clock() - st) / 1000.0);
		//printf("%d th node: %d, covered paths: %.3lf\n", k + 1, 	largestNodeID, largestIncreasedBC);
	}// end for (int k = 0; k < K; ++k)
	for (i = 0; i < n; ++i) // may be less than K nodes are chosen
	{
		if (foundNodes.size() >= K) break;
		if (isFoundforGBC[i] == true) continue;
		foundNodes.push_back(i); // arbitraries 
	}

	totalGBC = (totalGBC / L);
	for (i = 0; i < L; ++i) delete[]sampleShortestPathsArray[i];

	for (i = 0; i < n; ++i) {
		delete[]originalGraph[i];
		delete[] transposeGraph[i];
		delete[] predecessorsAdjLists[i];
		delete[] predecessorAdjListsBackwards[i];
	}
	delete[]originalGraphSize;
	delete[]transposeGraphSize;
	delete[]predecessorsAdjListsSize;
	delete[]predecessorAdjListsBackwardsSize;

	return totalGBC; // normalized GBC
}

double Graph::evaluateGBCbySampleShoretstPaths(vector<int>& foundNodes, vector<PAIR>& SamplePairs)
{
	clock_t  st = clock();
	int K = foundNodes.size();
	
	assert(K > 0 );
	int L = SamplePairs.size();
	assert(L <= numSamples);
	PAIR onePair;
	//printf("num of sampled paths: %d for evaluation\n", L);

	int i, j, k;
	int destID;
	bool debug = false;
	nonUniformArray<int>  predecessorsAdjLists; // record predecessors in the shortest paths
	vector<int> degrees; // degree of each node 
	degrees.reserve(n);
	for (i = 0; i < n; ++i) degrees.push_back(adjLists.size(i));
	predecessorsAdjLists.allocateMemory(degrees);

	int src, dst;
	int u, v;
	// the number of nodes that their shortest distances to src have been found
	//int numDistancesFound = 0;
	int sz;
	int* adj;

	double totalGBC = 0;

	//printf("num nodes: %d\n", n);
	int newDistance;
	double p;
	double delta;

	//int* sampleShortestPathsArray[numSamples];
	const int maxNodesInPath = 100;
	for (i = totalSampledShortestPathsSetT; i < L; ++i) sampleShortestPathsArraySetT[i] = new int[maxNodesInPath];

	degrees.clear();
	for (i = 0; i < n; ++i) degrees.push_back(0); // degrees in the transpose graph
	for (i = 0; i < n; ++i)
	{
		sz = adjLists.size(i);
		adj = &adjLists.access(i, 0);
		for (j = 0; j < sz; ++j)
			degrees[adj[j]]++;
	}
	nonUniformArray<int> transposeGraph;
	transposeGraph.allocateMemory(degrees);
	for (i = 0; i < n; ++i)
	{
		sz = adjLists.size(i);
		adj = &adjLists.access(i, 0);
		for (j = 0; j < sz; ++j)
			transposeGraph.push_back(adj[j], i);
	}
	nonUniformArray<int>  predecessorAdjListsBackwards;
	predecessorAdjListsBackwards.allocateMemory(degrees);
	//printf("original graph:\n"); adjLists.printArray();
	//printf("\n transpose graph:\n"); transposeGraph.printArray();


	clock_t ft = st;
	int middleNode;
	int tempNode;

	NodeInTwoQuques nodeU, nodeV;
	int stopBFSdistanceForward, stopBFSdistanceBackward;
	vector<EdgeWithWeight> touchEdges;
	double totalWeight;
	double sumWeight;

	vector<int> nodesInForwardBFS, nodesInBackwardBFS;
	nodesInForwardBFS.reserve(n);
	nodesInBackwardBFS.reserve(n);

	for (j = 0; j < n; ++j) predecessorsAdjLists.clear(j);
	for (j = 0; j < n; ++j) predecessorAdjListsBackwards.clear(j);
	for (j = 0; j < n; ++j) visited[j] = false;
	for (j = 0; j < n; ++j) visitedBackwards[j] = false;
	for (j = 0; j < n; ++j) numShortestPathsBefore[j] = 0;
	for (j = 0; j < n; ++j) numShortestPathsBackward[j] = 0;

	for (i = totalSampledShortestPathsSetT; i <  L; ++i) // randomly generate L shortest paths
	{
		//if (i == 40) debug = true;
		//else debug = false;
		//debug = true;

		onePair = SamplePairs[i];
		src = onePair.s;
		dst = onePair.t;

		/*if ((i + 1) % 100000 == 0) {
			printf("L: %d, find the %d th sample, time used: %.2lf, delta time: %.2lf\n", L, i + 1, (clock() - st) / 1000.0,
				(clock() - ft) / 1000.0);
			ft = clock();
		}*/
		

		// find the shortest paths from src 
		sz = nodesInForwardBFS.size();
		for (j = 0; j < sz; ++j) {
			v = nodesInForwardBFS[j];
			predecessorsAdjLists.clear(v);
			visited[v] = false;
			numShortestPathsBefore[v] = 0;
		}
		visited[src] = true;
		numShortestPathsBefore[src] = 1;
		nodesInForwardBFS.clear();
		nodesInForwardBFS.push_back(src);

		// find the shortest paths from dst 
		sz = nodesInBackwardBFS.size();
		for (j = 0; j < sz; ++j) {
			v = nodesInBackwardBFS[j];
			predecessorAdjListsBackwards.clear(v);
			visitedBackwards[v] = false;
			numShortestPathsBackward[v] = 0;
		}
		visitedBackwards[dst] = true;
		numShortestPathsBackward[dst] = 1;
		nodesInBackwardBFS.clear();
		nodesInBackwardBFS.push_back(dst);

		distances[src] = 0;
		distancesBackwards[dst] = 0;
		biQueue_clear(); // clear the queue
		biQueue_push(NodeInTwoQuques(src, 1));  //first queue
		biQueue_push(NodeInTwoQuques(dst, 2));  // second queue
		stopBFSdistanceForward = stopBFSdistanceBackward = n - 1;
		touchEdges.clear();
		//if (debug) printf("\npair(%d, %d)******************\n", src, dst);

		while (biQueue_size() > 0) // bfs
		{
			nodeU = biQueue_pop();
			u = nodeU.v;
			//if (debug) printf("exploring node %d ", u);

			if (1 == nodeU.qID) // forward search
			{
				//if (debug) printf("forward \n");
				if (distances[u] >= stopBFSdistanceForward) continue;

				sz = adjLists.size(u);
				adj = &adjLists.access(u, 0);
				newDistance = distances[u] + 1;

				for (j = 0; j < sz; ++j)
				{
					v = adj[j];
					if (false == visited[v])
					{
						visited[v] = true;
						distances[v] = newDistance;
						biQueue_push(NodeInTwoQuques(v, 1));
						nodesInForwardBFS.push_back(v);
						//if (debug) printf("node %d is added to queue\n", v);
					}
					if (distances[v] == newDistance) {
						predecessorsAdjLists.push_back(v, u);
						numShortestPathsBefore[v] += numShortestPathsBefore[u];
					}

					if (true == visitedBackwards[v])
					{
						stopBFSdistanceForward = newDistance;
						stopBFSdistanceBackward = distancesBackwards[v];
						touchEdges.push_back(EdgeWithWeight(u, v, 1, 0));  // forward edge
						//if (debug) printf("at node %d, found a middle node: %d\n", u, v);
					}
				}
			}
			else { // backward search
				//if (debug) printf("backward \n");
				if (distancesBackwards[u] >= stopBFSdistanceBackward) continue;
				sz = transposeGraph.size(u);
				adj = &transposeGraph.access(u, 0);
				newDistance = distancesBackwards[u] + 1;
				for (j = 0; j < sz; ++j)
				{
					v = adj[j];
					if (false == visitedBackwards[v])
					{
						visitedBackwards[v] = true;
						distancesBackwards[v] = newDistance;
						biQueue_push(NodeInTwoQuques(v, 2));
						nodesInBackwardBFS.push_back(v);
						//if (debug) printf("node %d is added to queue\n", v);
					}
					if (distancesBackwards[v] == newDistance) {
						predecessorAdjListsBackwards.push_back(v, u);
						numShortestPathsBackward[v] += numShortestPathsBackward[u];
					}

					if (true == visited[v])
					{
						stopBFSdistanceForward = distances[v];
						stopBFSdistanceBackward = newDistance;
						touchEdges.push_back(EdgeWithWeight(u, v, 2, 0)); // backward edge
						//if (debug) printf("at node %d, found a middle node: %d\n", u, v);
					}
				}
			}
		}
		destID = onePair.t;
		nodesInEachSampleSetT[i] = 0;
		if (touchEdges.size() == 0) {
			//printf("pair: %d, %d, not reachable\n", src, destID);
			//abort();
			//nodesInEachSample[i] = 0;
			continue;
		}

		totalWeight = 0;
		for (j = 0; j < touchEdges.size(); ++j)
		{
			if (1 == touchEdges[j].direction) // forward edge
			{
				u = touchEdges[j].u;
				v = touchEdges[j].v;
			}
			else { // backward edge
				u = touchEdges[j].v;
				v = touchEdges[j].u;
			}
			touchEdges[j].weight = double(numShortestPathsBefore[u]) * numShortestPathsBackward[v];
			totalWeight += touchEdges[j].weight;
			//if (debug) printf("touch edge: (%d, %d), weight: %.1lf\n", u, v, touchEdges[j].weight);
		}
		p = Uniform();
		//if (debug) printf("p: %.3lf\n", p);
		p *= totalWeight;
		sumWeight = 0;

		for (j = 0; j < touchEdges.size(); ++j)
		{
			sumWeight += touchEdges[j].weight;
			if (sumWeight >= p) break;
		}
		assert(j != touchEdges.size());
		middleNode = touchEdges[j].v;
		//if (debug) printf("chosen edge: (% d, % d), weight: % .1lf\n", touchEdges[j].u, touchEdges[j].v, touchEdges[j].weight);

		/*if (debug)
		{
			printf("pair: (%d, %d)\n predecessorsAdjLists:\n", src, onePair.t);
			predecessorsAdjLists.printArray();
			printf("predecessorAdjListsBackwords:\n");
			predecessorAdjListsBackwards.printArray();
			printf("middle node: %d\n", middleNode);
		}*/

		tempNode = middleNode;
		while (tempNode != src)
		{
			//assert(tempNode >= 0 && tempNode < n);
			//if (nodesInEachSample[i] >= maxNodesInPath) abort();
			sampleShortestPathsArraySetT[i][nodesInEachSampleSetT[i]] = tempNode;
			nodesInEachSampleSetT[i]++;
			p = Uniform();
			sz = predecessorsAdjLists.size(tempNode);
			/*if (sz == 0) {
				printf("error: pair: %d, %d\npredecessorsAdjLists:\n", onePair.s, onePair.t);
				//predecessorsAdjLists.printArray();
				printf("%d th sample\n", i + 1);
				abort();
			}*/
			delta = 1.0 / sz;
			j = int(floor(p / delta)); // round bin gaming
			if (j == sz) j--;
			//assert(j >= 0 && j < sz);
			//if (debug) printf("---INSERT node %d\n", tempNode);

			tempNode = predecessorsAdjLists.access(tempNode, j);
		}
		//if (nodesInEachSample[i] >= maxNodesInPath) abort();
		sampleShortestPathsArraySetT[i][nodesInEachSampleSetT[i]] = src;
		nodesInEachSampleSetT[i]++;
		//if (debug) printf("---INSERT node %d\n", src);

		tempNode = middleNode;
		if (tempNode == dst) continue;
		p = Uniform();
		sz = predecessorAdjListsBackwards.size(tempNode);
		delta = 1.0 / sz;
		j = int(floor(p / delta)); // round bin gaming
		if (j == sz) j--;
		//assert(j >= 0 && j < sz);
		tempNode = predecessorAdjListsBackwards.access(tempNode, j);

		while (tempNode != dst)
		{
			//assert(tempNode >= 0 && tempNode < n);
			//if (nodesInEachSample[i] >= maxNodesInPath) abort();
			sampleShortestPathsArraySetT[i][nodesInEachSampleSetT[i]] = tempNode;
			nodesInEachSampleSetT[i]++;
			p = Uniform();
			sz = predecessorAdjListsBackwards.size(tempNode);
			/*if (sz == 0) {
				printf("error: pair: %d, %d \n predecessorAdjListsBackwards:\n", onePair.s, onePair.t);
				//predecessorAdjListsBackwards.printArray();
				printf("%d th sample\n", i + 1);
				abort();
			}*/
			delta = 1.0 / sz;
			j = int(floor(p / delta)); // round bin gaming
			if (j == sz) j--;
			//if (j < 0) j = 0;
			//assert(j >= 0 && j < sz);
			//if (debug) printf("---INSERT node %d\n", tempNode);

			tempNode = predecessorAdjListsBackwards.access(tempNode, j);
		}
		//if (nodesInEachSample[i] >= maxNodesInPath) abort();
		sampleShortestPathsArraySetT[i][nodesInEachSampleSetT[i]] = dst;
		nodesInEachSampleSetT[i]++;
		//if (debug) printf("---INSERT node %d\n", dst);
		/*if (debug)
		{
			printf("%d th sample: ", i + 1);
			for (j = 0; j < nodesInEachSample[i]; ++j)
				printf("%d,", sampleShortestPathsArray[i][j]);
			printf(", middle Point: %d\n", middleNode);
		}*/
	}
	totalSampledShortestPathsSetT = L;
	
	//printf("bi BFS numEdegesExplord: %.0lf\n", numEdegesExplord);

	for (i = 0; i < n; ++i) isFoundforGBC[i] = false;
	for (i = 0; i < K; ++i) isFoundforGBC[foundNodes[i]] = true;

	bool isCovered;
	double numCoveredPaths = 0;
	for (i = 0; i < L; ++i)
	{
		isCovered = false;
		for (j = 0; j < nodesInEachSampleSetT[i]; ++j)
		{
			v = sampleShortestPathsArraySetT[i][j];
			if (true == isFoundforGBC[v]) {
				isCovered = true;
				break;
			}
		}
		if (true == isCovered)
			numCoveredPaths++;
	}
	totalGBC = numCoveredPaths / L; 
	
	//for (i = 0; i < L; ++i) delete[]sampleShortestPathsArray[i];
	return totalGBC; // normalized GBC
}

void Graph::generateRandomPairs(int n, int L, vector<PAIR>& randomizedPairs)
{// Post: generate L pairs, with indices between 0 and n-1
	assert(n > 0);
	randomizedPairs.clear();
	randomizedPairs.reserve(L);

	int i, j;
	int s, t;
	double p;
	PAIR pair;
	for (i = 0; i < L; ++i)
	{
		p = Uniform();
		s = int(floor(p * n));
		if (s == n) s--;
		pair.s = s;
		while (true) {
			p = Uniform();
			t = int(floor(p * n));
			if (t == n) t--;
			if (t != s) break;
		}
		pair.t = t;
		randomizedPairs.push_back(pair);
	}
	/*for (i = 0; i < L; ++i)
	{
		s = randomizedPairs[i].s;
		t = randomizedPairs[i].t;
		assert(s >= 0 && s < n);
		assert(t >= 0 && t < n && s!= t);
	}*/
}


double Graph::findTopKGBCAdaptiveSampling(int K, double epsilon, double errorProb, vector<int>& foundNodes, int& sampledPaths)
{
	assert(K > 0 && K < n);
	assert(epsilon > 0 && epsilon < 0.632);
	assert(errorProb > 0 && errorProb = 1);
	foundNodes.clear();
	foundNodes.reserve(K);

	double e = 2.71828;
	double alpha = epsilon / (2.0 - 1 / e);
	double c2 = (2 + alpha) / (alpha * alpha);
	double b = 3 * c2 + 2 + sqrt(18*c2+4);
	b /= 3 * c2 - 2;
	if (b < 1.1) b = 1.1;
	//b = 2; // for testing

	double Qmax = log(double(n) * (n - 1.0));
	Qmax = ceil( Qmax / log(b) );
	//printf("\nn: %d, base b: %.3lf, Qmax: %.0lf\n", n, b, Qmax);

	double d1;
	 d1 = log(2.0 / errorProb) + log(Qmax);

	double cnt = 0;
	int q;
	
	double Lq;
	double gq = 1.0;
	double theta = d1 * c2; // (ln 2/gamma+ln Qmax ) (2+alpha)/alpha^2
	Lq = theta;
	vector<PAIR> randomPairs;
	vector<int> curFoundNodes;
	vector<double> foundTime;
	double curBiasedGBC, curUnbiasedGBC;
	double totalSamples = 0;

	double epsilon1,  epsilonSum;
	double beta;
	double betaMax = epsilon / (1.0 - 1.0 / e);
	double timesLessThanOPT;
	double c1;

	totalSampledShortestPathsSetS = 0;
	totalSampledShortestPathsSetT = 0;
	for (q = 1; q <= Qmax; ++q)
	{
		gq /= b;
		Lq *= b;

		generateRandomPairs(n, int(Lq), randomPairs);
		curBiasedGBC = rankByGroupBCWithSampleShortestPathsBiBFS(randomPairs, curFoundNodes, K, foundTime);

		generateRandomPairs(n, int(Lq), randomPairs);
		curUnbiasedGBC = evaluateGBCbySampleShoretstPaths(curFoundNodes, randomPairs);
		beta = 1.0 - curUnbiasedGBC / curBiasedGBC;

		//printf("%dth try, biased gbc: %.4lf, unbiased gbc: %.4lf, gq: %.4lf, L: %d, totalSampledShortestPaths: %d\n", q, curBiasedGBC, curUnbiasedGBC, gq, (int)Lq, totalSampledShortestPathsSetS);
		if (curUnbiasedGBC >= gq) cnt++;
		if (cnt >= 2) {
			timesLessThanOPT = pow(b, cnt - 2);
			c1 = log(4.0 / errorProb) / (theta * timesLessThanOPT);
			epsilon1 = (2*c1/3.0 + sqrt(4*c1 * c1/9.0 + 8 * c1)) / 2.0;
			//epsilon2 = sqrt(2*log(4.0/errorProb)/(theta * timesLessThanOPT) );
			epsilonSum = beta * (1 - 1.0 / e) * (1 - epsilon1) + (2 - 1.0/e) * epsilon1;
			betaMax = 1.0 - (1 - 1 / e - epsilon + epsilon1) / ((1 - 1 / e) * (1 - epsilon1));
			//printf("cnt: %.0lf, epsilon1: %.4lf,  beta: %.4lf, betaMax: %4lf, epsilonSum: %.4lf\n", cnt,	epsilon1,  beta, betaMax, epsilonSum);
			if (epsilonSum > epsilon) continue;

			for (int k = 0; k < totalSampledShortestPathsSetS; ++k)
				delete[]sampleShortestPathsArraySetS[k];
			for (int k = 0; k < totalSampledShortestPathsSetT; ++k)
				delete[]sampleShortestPathsArraySetT[k];
			sampledPaths = totalSampledShortestPathsSetS + totalSampledShortestPathsSetT;
			for (int k = 0; k < K; ++k) foundNodes.push_back(curFoundNodes[k]);
			//sampledPaths = totalSamples;
			
			return curUnbiasedGBC;
		}
	}
}


/*
double Graph::findTopKGBCbyKDD23(int K, double epsilon, double errorProb, vector<int>& foundNodes, int& sampledPaths)
{
	assert(K > 0 && K < n);
	assert(epsilon > 0 && epsilon < 0.632);
	assert(errorProb > 0 && errorProb = 1);
	foundNodes.clear();
	foundNodes.reserve(K);

	double e = 2.71828;
	double alpha = epsilon / (2.0 - 1 / e);
	double c1 = (2 + alpha) / (alpha * alpha);

	double Qmax = ceil(log(double(n) * (n - 1.0)));
	double c2 = log(2.0 / errorProb) + log(Qmax);

	double cnt = 0;
	int q;
	double Lq;
	double gq = 1.0;
	Lq = (c2 + K * log( 3.0 * K)) * c1; // (ln 1/gamma+ln Qmax + K ln n) (2+alpha)/alpha^2
	vector<PAIR> randomPairs;
	vector<int> curFoundNodes;
	vector<double> foundTime;
	double curGBC;
	double totalSamples = 0;
	for (q = 1; q <= Qmax; ++q)
	{
		gq /= 2.0;
		Lq *= 2.0;
		totalSamples += Lq;

		generateRandomPairs(n, int(Lq), randomPairs);
		curGBC = rankByGroupBCWithSampleShortestPathsBiBFSunbiased(randomPairs, curFoundNodes, K, foundTime);

		//printf("%dth try, gbc: %.4lf, gq: %.4lf, Lq: %d\n", q, curGBC, gq, (int)Lq);
		if (curGBC >= gq) cnt++;
		if (cnt >= 2) {
			for (int k = 0; k < K; ++k) foundNodes.push_back(curFoundNodes[k]);
			sampledPaths = totalSamples;
			return curGBC;
		}
	}
}
/**/

double Graph::findTopKGBCbyKDD16(int K, double epsilon, double errorProb, vector<int>& foundNodes, int &sampledPaths)
{
	assert(K > 0 && K < n);
	assert(epsilon > 0 && epsilon < 0.632);
	assert(errorProb > 0 && errorProb = 1);
	foundNodes.clear();
	foundNodes.reserve(K);

	double e = 2.71828;
	double alpha = epsilon / (2.0 - 1 / e);
	double c1 = (2 + alpha) / (alpha * alpha);

	double Qmax = ceil(log(double(n) * (n - 1.0)));
	double c2 = log(2.0 / errorProb) + log(Qmax);

	double cnt = 0;
	int q;
	double Lq;
	double gq = 1.0;
	Lq = (c2 + K * log(double(n))) * c1; // (ln 1/gamma+ln Qmax + K ln n) (2+alpha)/alpha^2
	vector<PAIR> randomPairs;
	vector<int> curFoundNodes;
	vector<double> foundTime;
	double curGBC;
	double totalSamples = 0;

	for (q = 1; q <= Qmax; ++q)
	{
		gq /= 2.0;
		Lq *= 2.0;
		totalSamples += Lq;

		generateRandomPairs(n, int(Lq), randomPairs);
		curGBC = rankByGroupBCWithSampleShortestPathsBiBFS(randomPairs, curFoundNodes, K, foundTime);
		
		//printf("%dth try, gbc: %.4lf, gq: %.4lf, Lq: %d\n", q, curGBC, gq, (int)Lq);
		if (curGBC >= gq) cnt++;
		if (cnt >= 2){
			for (int k = 0; k < K; ++k) foundNodes.push_back(curFoundNodes[k]);
			sampledPaths = totalSamples;
			return curGBC;
		}
	}
}



double Graph::rankByGroupBCWithSampleShortestPathsBiBFS(vector<PAIR>& randomPairs, vector<int>& foundNodes, int K, vector<double>& foundTimes)
{ // find the sample with bidirectional BFS
	clock_t  st = clock();

	foundNodes.clear();
	foundNodes.reserve(K);
	foundTimes.clear();
	assert(K > 0 && K <= n);
	int L = randomPairs.size();
	assert(L <= numSamples);
	PAIR onePair;

	int i, j, k;
	int destID;
	bool debug = false;
	nonUniformArray<int>  predecessorsAdjLists; // record predecessors in the shortest paths
	vector<int> degrees; // degree of each node 
	degrees.reserve(n);
	for (i = 0; i < n; ++i) degrees.push_back(adjLists.size(i));
	predecessorsAdjLists.allocateMemory(degrees);

	int src, dst;
	int u, v;
	// the number of nodes that their shortest distances to src have been found
	//int numDistancesFound = 0;
	int sz;
	int* adj;

	double totalGBC = 0;
	double largestIncreasedBC;
	int largestNodeID;

	//printf("num nodes: %d\n", n);
	int newDistance;
	double p;
	double delta;

	
	const int maxNodesInPath = 100;
	for (i = totalSampledShortestPathsSetS; i < L; ++i) sampleShortestPathsArraySetS[i] = new int[maxNodesInPath];

	degrees.clear();
	for (i = 0; i < n; ++i) degrees.push_back(0); // degrees in the transpose graph
	for (i = 0; i < n; ++i)
	{
		sz = adjLists.size(i);
		adj = &adjLists.access(i, 0);
		for (j = 0; j < sz; ++j)
			degrees[adj[j]]++;
	}
	nonUniformArray<int> transposeGraph;
	transposeGraph.allocateMemory(degrees);
	for (i = 0; i < n; ++i)
	{
		sz = adjLists.size(i);
		adj = &adjLists.access(i, 0);
		for (j = 0; j < sz; ++j)
			transposeGraph.push_back(adj[j], i);
	}
	nonUniformArray<int>  predecessorAdjListsBackwards;
	predecessorAdjListsBackwards.allocateMemory(degrees);
	//printf("original graph:\n"); adjLists.printArray();
	//printf("\n transpose graph:\n"); transposeGraph.printArray();

	int middleNode;
	int tempNode;

	NodeInTwoQuques nodeU, nodeV;
	int stopBFSdistanceForward, stopBFSdistanceBackward;
	vector<EdgeWithWeight> touchEdges;
	touchEdges.reserve(100);
	double totalWeight;
	double sumWeight;

	clock_t ft = st;
	//printf("time used for initialization: %.3lf\n", (ft - st) / 1000.0);

	vector<int> nodesInForwardBFS, nodesInBackwardBFS;
	nodesInForwardBFS.reserve(n); 
	nodesInBackwardBFS.reserve(n);

	for (j = 0; j < n; ++j) predecessorsAdjLists.clear(j);
	for (j = 0; j < n; ++j) predecessorAdjListsBackwards.clear(j);
	for (j = 0; j < n; ++j) visited[j] = false;
	for (j = 0; j < n; ++j) visitedBackwards[j] = false;
	for (j = 0; j < n; ++j) numShortestPathsBefore[j] = 0;
	for (j = 0; j < n; ++j) numShortestPathsBackward[j] = 0;

	for (i = totalSampledShortestPathsSetS; i < L; ++i) // randomly generate L shortest paths
	{
		//if (i == 40) debug = true;
		//else debug = false;
		//debug = true;

		onePair = randomPairs[i];
		src = onePair.s;
		dst = onePair.t;

		if ((i+1) % 1000000 == 0)printf("L: %d, find the %d th sample, time used: %.2lf\n", L, i + 1, (clock() - st) / 1000.0);


		// find the shortest paths from src 
		sz = nodesInForwardBFS.size();
		for (j = 0; j < sz; ++j) {
			v = nodesInForwardBFS[j];
			predecessorsAdjLists.clear(v);
			visited[v] = false;
			numShortestPathsBefore[ v ] = 0;
		}
		visited[src] = true;
		numShortestPathsBefore[src] = 1;
		nodesInForwardBFS.clear();
		nodesInForwardBFS.push_back(src);
		
		// find the shortest paths from dst 
		sz = nodesInBackwardBFS.size();
		for (j = 0; j < sz; ++j) {
			v = nodesInBackwardBFS[j];
			predecessorAdjListsBackwards.clear(v);
			visitedBackwards[v] = false;
			numShortestPathsBackward[v] = 0;
		}
		visitedBackwards[dst] = true;
		numShortestPathsBackward[dst] = 1;
		nodesInBackwardBFS.clear();
		nodesInBackwardBFS.push_back(dst);

		distances[src] = 0;
		distancesBackwards[dst] = 0;
		biQueue_clear(); // clear the queue
		biQueue_push(NodeInTwoQuques(src, 1));  //first queue
		biQueue_push(NodeInTwoQuques(dst, 2));  // second queue
		stopBFSdistanceForward = stopBFSdistanceBackward = n - 1;
		touchEdges.clear();
		if( debug ) printf("\npair(%d, %d)******************\n", src, dst);

		while (biQueue_size() > 0) // bfs
		{
			nodeU = biQueue_pop();
			u = nodeU.v;
			//if (debug) printf("exploring node %d ", u);

			if (1 == nodeU.qID ) // forward search
			{
				//if (debug) printf("forward \n");
				if (distances[u] >= stopBFSdistanceForward) continue;

				sz = adjLists.size(u);
				adj = adjLists.startAddress(u);
				newDistance = distances[u] + 1;
				for (j = 0; j < sz; ++j)
				{
					v = adj[j];
					if (false == visited[v])
					{
						visited[v] = true;
						distances[v] = newDistance;
						biQueue_push(NodeInTwoQuques(v,1));
						nodesInForwardBFS.push_back(v);
						//if (debug) printf("node %d is added to queue\n", v);
					}
					if (distances[v] == newDistance) {
						predecessorsAdjLists.push_back(v, u);
						numShortestPathsBefore[v] += numShortestPathsBefore[u];
					}
						
					if (true == visitedBackwards[v])
					{
						stopBFSdistanceForward = newDistance;
						stopBFSdistanceBackward = distancesBackwards[v];
						touchEdges.push_back(EdgeWithWeight(u, v, 1, 0));  // forward edge
						//if (debug) printf("at node %d, found a middle node: %d\n", u, v);
					}
				}
			}
			else { // backward search
				//if (debug) printf("backward \n");
				if (distancesBackwards[ u ] >= stopBFSdistanceBackward) continue;
				sz = transposeGraph.size(u);
				adj = transposeGraph.startAddress(u);
				newDistance = distancesBackwards[u] + 1;
				for (j = 0; j < sz; ++j)
				{
					v = adj[j];
					if (false == visitedBackwards[v])
					{
						visitedBackwards[v] = true;
						distancesBackwards[v] = newDistance;
						biQueue_push(NodeInTwoQuques(v,2));
						nodesInBackwardBFS.push_back(v);
						//if (debug) printf("node %d is added to queue\n", v);
					}
					if (distancesBackwards[v] == newDistance) {
						predecessorAdjListsBackwards.push_back(v, u);
						numShortestPathsBackward[v] += numShortestPathsBackward[u];
					}
						
					if (true == visited[v])
					{
						stopBFSdistanceForward = distances[v];
						stopBFSdistanceBackward = newDistance;
						touchEdges.push_back(EdgeWithWeight(u, v, 2, 0)); // backward edge
						//if (debug) printf("at node %d, found a middle node: %d\n", u, v);	
					}
				}
			}
		}
		destID = onePair.t;
		nodesInEachSampleSetS[i] = 0;
		if (touchEdges.size() == 0) {
			//printf("pair: %d, %d, not reachable\n", src, destID);
			//abort();
			continue;
		}

		totalWeight = 0;
		for (j = 0; j < touchEdges.size(); ++j)
		{
			if (1 == touchEdges[j].direction) // forward edge
			{
				u = touchEdges[j].u;
				v = touchEdges[j].v;
			}else { // backward edge
				u = touchEdges[j].v;
				v = touchEdges[j].u;
			}
			touchEdges[j].weight = double(numShortestPathsBefore[u]) * numShortestPathsBackward[v];
			totalWeight += touchEdges[j].weight;
			//if (debug) printf("touch edge: (%d, %d), weight: %.1lf\n", u, v, touchEdges[j].weight);
		}
		p = Uniform();
		//if (debug) printf("p: %.3lf\n", p);
		p *= totalWeight;
		sumWeight = 0;
		for (j = 0; j < touchEdges.size(); ++j)
		{
			sumWeight += touchEdges[j].weight;
			if (sumWeight >= p) break;
		}
		//assert(j != touchEdges.size());
		middleNode = touchEdges[j].v;
		//if (debug) printf("chosen edge: (% d, % d), weight: % .1lf\n", touchEdges[j].u, touchEdges[j].v, touchEdges[j].weight);

		/*if (debug)
		{
			printf("pair: (%d, %d)\n predecessorsAdjLists:\n", src, onePair.t);
			predecessorsAdjLists.printArray();
			printf("predecessorAdjListsBackwords:\n");
			predecessorAdjListsBackwards.printArray();
			printf("middle node: %d\n", middleNode);
		}*/

		tempNode = middleNode;
		while (tempNode != src)
		{
			//assert(tempNode >= 0 && tempNode < n);
			//if (nodesInEachSample[i] >= maxNodesInPath) abort();
			sampleShortestPathsArraySetS[i][nodesInEachSampleSetS[i]] = tempNode;
			nodesInEachSampleSetS[i]++;
			p = Uniform();
			sz = predecessorsAdjLists.size(tempNode);
			/*if (sz == 0) {
				printf("error: pair: %d, %d\npredecessorsAdjLists:\n", onePair.s, onePair.t);
				//predecessorsAdjLists.printArray();
				printf("%d th sample\n", i + 1);
				abort();
			}*/
			delta = 1.0 / sz;
			j = int(floor(p / delta)); // round bin gaming
			if (j == sz) j--;
			//assert(j >= 0 && j < sz);
			//if (debug) printf("---INSERT node %d\n", tempNode);

			tempNode = predecessorsAdjLists.access(tempNode, j);
		}
		//if (nodesInEachSample[i] >= maxNodesInPath) abort();
		sampleShortestPathsArraySetS[i][nodesInEachSampleSetS[i]] = src;
		nodesInEachSampleSetS[i]++;
		//if (debug) printf("---INSERT node %d\n", src);

		tempNode = middleNode;
		if (tempNode == dst) continue;
		p = Uniform();
		sz = predecessorAdjListsBackwards.size(tempNode);
		delta = 1.0 / sz;
		j = int(floor(p / delta)); // round bin gaming
		if (j == sz) j--;
		//assert(j >= 0 && j < sz);
		tempNode = predecessorAdjListsBackwards.access(tempNode, j);

		while (tempNode != dst)
		{
			//assert(tempNode >= 0 && tempNode < n);
			//if (nodesInEachSample[i] >= maxNodesInPath) abort();
			sampleShortestPathsArraySetS[i][nodesInEachSampleSetS[i]] = tempNode;
			nodesInEachSampleSetS[i]++;
			p = Uniform();
			sz = predecessorAdjListsBackwards.size(tempNode);
			/*if (sz == 0) {
				printf("error: pair: %d, %d \n predecessorAdjListsBackwards:\n", onePair.s, onePair.t);
				//predecessorAdjListsBackwards.printArray();
				printf("%d th sample\n", i + 1);
				abort();
			}*/
			delta = 1.0 / sz;
			j = int(floor(p / delta)); // round bin gaming
			if (j == sz) j--;
			//if (j < 0) j = 0;
			//assert(j >= 0 && j < sz);
			//if (debug) printf("---INSERT node %d\n", tempNode);

			tempNode = predecessorAdjListsBackwards.access(tempNode, j);
		}
		//if (nodesInEachSample[i] >= maxNodesInPath) abort();
		sampleShortestPathsArraySetS[i][nodesInEachSampleSetS[i]] = dst;
		nodesInEachSampleSetS[i]++;
		//if (debug) printf("---INSERT node %d\n", dst);
		/*if (debug)
		{
			printf("%d th sample: ", i + 1);
			for (j = 0; j < nodesInEachSample[i]; ++j)
				printf("%d,", sampleShortestPathsArray[i][j]);
			printf(", middle Point: %d\n", middleNode);
		}*/
	}
	foundTimes.push_back((clock() - st) / 1000.0);
	totalSampledShortestPathsSetS = L;

	for (i = 0; i < n; ++i) isFoundforGBC[i] = false;
	int* pointerOnePath;
	for (k = 0; k < K; ++k)
	{
		//if ( (k+1) % 10 == 0) printf("finding the %d node...\n", k + 1);
		for (i = 0; i < n; ++i) betweenCentralities[i] = 0;
		for (i = 0; i < L; ++i)
		{
			pointerOnePath = sampleShortestPathsArraySetS[i];
			for (j = 0; j < nodesInEachSampleSetS[i]; ++j)
				if (isFoundforGBC[pointerOnePath[j]] == true) break;
			if (j != nodesInEachSampleSetS[i]) continue; // already covered
			for (j = 0; j < nodesInEachSampleSetS[i]; ++j)
				betweenCentralities[pointerOnePath[j]]++;
		}
		largestIncreasedBC = 0;
		for (i = 0; i < n; ++i)
		{
			//if (betweenCentralities[i] > 0)
				//assert(isFoundforGBC[i] == false);
			if (betweenCentralities[i] > largestIncreasedBC)
			{
				largestIncreasedBC = betweenCentralities[i];
				largestNodeID = i;
			}
		}
		if (largestIncreasedBC < 0.5) break;// no nodes are found
		isFoundforGBC[largestNodeID] = true;
		foundNodes.push_back(largestNodeID);
		totalGBC += largestIncreasedBC;

		if ((k + 1) % 10 == 0) foundTimes.push_back((clock() - st) / 1000.0);
		//printf("%d th node: %d, covered paths: %.3lf\n", k + 1, 	largestNodeID, largestIncreasedBC);
	}// end for (int k = 0; k < K; ++k)
	for (i = 0; i < n; ++i) // may be less than K nodes are chosen
	{
		if (foundNodes.size() >= K) break;
		if (isFoundforGBC[i] == true) continue;
		foundNodes.push_back(i); // arbitrarily
	}
	//printf("time used for finding: %.3lf\n", (clock() - ft) / 1000.0);

	//for (i = 0; i < L; ++i) delete[]sampleShortestPathsArray[i];

	totalGBC /= L;
	return totalGBC; // normalized GBC
}


/*
double Graph::rankByGroupBCWithSampleShortestPaths(vector<PAIR>& randomPairs, vector<int>& foundNodes, int K, vector<double>& foundTimes)
{
	foundNodes.clear();
	foundNodes.reserve(K);
	foundTimes.clear();
	assert(K > 0 && K <= n);
	int L = randomPairs.size();
	assert(L <= numSamples  );
	PAIR onePair;

	int i, j, k;
	int destID;
	bool debug = false;
	nonUniformArray<int>  predecessorsAdjLists; // record predecessors in the shortest paths
	vector<int> degrees; // degree of each node 
	degrees.reserve(n);
	for (i = 0; i < n; ++i) degrees.push_back(adjLists.size(i));
	predecessorsAdjLists.allocateMemory(degrees);
	int src;
	int u, v;
	// the number of nodes that their shortest distances to src have been found
	//int numDistancesFound = 0;
	int sz;
	int* adj;

	double totalGBC = 0;
	double largestIncreasedBC;
	int largestNodeID;

	//printf("num nodes: %d\n", n);
	int newDistance;
	double p;
	double delta;

	int* sampleShortestPathsArray[numSamples];
	const int maxNodesInPath = 100;
	for (i = 0; i < L; ++i) sampleShortestPathsArray[i] = new int[maxNodesInPath];

	clock_t  st = clock();
	clock_t ft = st;
	double numEdegesExplord = 0;
	for (i = 0; i < L; ++i) // randomly generate L shortest paths
	{
		onePair = randomPairs[i];
		src = onePair.s;
		
		if (i % 1000 == 0) {
			printf("L: %d, find the %d th sample, time used: %.2lf, delta time: %.2lf\n", L, i + 1, (clock()-st)/1000.0, 
				(clock() - ft) / 1000.0);
			ft = clock();
		}
		for (j = 0; j < n; ++j) predecessorsAdjLists.clear(j);
	
		// find the shortest paths from src 
		for (j = 0; j < n; ++j) visited[j] = false;
		visited[src] = true;
		distances[src] = 0;
		queue_clear(); // clear the queue
		queue_push(src);
		while (queue_size() > 0) // bfs
		{
			u = queue_pop();
			if (u == onePair.t) break;

			sz = adjLists.size(u);
			adj = &adjLists.access(u, 0);
			newDistance = distances[u] + 1;
			for (j = 0; j < sz; ++j)
			{
				v = adj[j];
				if (false == visited[v])
				{
					visited[v] = true;
					distances[v] = newDistance;
					queue_push(v);
				}
				if (distances[v] == newDistance) 
					predecessorsAdjLists.push_back(v, u);
			}
			numEdegesExplord += sz;
		}
		if (debug)
		{
			printf("pair: (%d, %d)\n predecessorsAdjLists:\n", src, onePair.t);
			predecessorsAdjLists.printArray();
		}
		destID = onePair.t;
		nodesInEachSample[i] = 0;
		if (false == visited[destID]){
			//printf("pair: %d, %d, not reachable\n", src, destID);
			//abort();
			sampleShortestPathsArray[i][0] = src;
			nodesInEachSample[i] = 1;
			continue;
 		}
		
		while (destID != src)
		{
			assert(destID >= 0 && destID < n);
			if (nodesInEachSample[i] >= maxNodesInPath) abort();
			sampleShortestPathsArray[i][ nodesInEachSample[i] ] = destID;
			nodesInEachSample[i]++;
			p = Uniform();
			sz = predecessorsAdjLists.size(destID);
			if (sz ==0) {
				//printf("error: pair: %d, %d\npredecessorsAdjLists:\n", onePair.s, onePair.t);
				//predecessorsAdjLists.printArray();
				abort();
			}
			delta = 1.0 / sz ;
			j = int (floor(p / delta)); // round bin gaming
			if (j == sz ) j--;
			//if (j < 0) j = 0;
			assert(j >= 0 && j < sz);

			destID = predecessorsAdjLists.access(destID, j);
		}
		if (nodesInEachSample[i] >= maxNodesInPath) abort();
		sampleShortestPathsArray[i][nodesInEachSample[i]] = src;
		nodesInEachSample[i]++;
		if (debug)
		{
			printf("%d th sample: ", i + 1);
			for (j = 0; j < nodesInEachSample[i]; ++j)
				printf("%d,", sampleShortestPathsArray[i][j]);
			printf("\n");
		}
	}
	foundTimes.push_back( (clock()-st)/1000.0);
	numEdegesExplord /= L;
	//printf("numEdegesExplord: %.0lf\n", numEdegesExplord);

	for (i = 0; i < n; ++i) isFoundforGBC[i] = false;
	for (k = 0; k < K; ++k)
	{
		//if( k%10 == 0) printf("finding the %d node...\n", k + 1);
		for (i = 0; i < n; ++i) betweenCentralities[i] = 0;
		for (i = 0; i < L; ++i)
		{
			for (j = 0; j < nodesInEachSample[i]; ++j)
				if (isFoundforGBC[ sampleShortestPathsArray[i][j] ] == true) break;
			if (j != nodesInEachSample[i]) continue; // already covered
			for (j = 0; j < nodesInEachSample[i]; ++j)
				betweenCentralities[ sampleShortestPathsArray[i][j] ]++;
		}
		largestIncreasedBC = 0;
		for (i = 0; i < n; ++i)
		{
			if (betweenCentralities[i] > 0)
				assert(isFoundforGBC[i] == false);
			if (betweenCentralities[i] > largestIncreasedBC)
			{
				largestIncreasedBC = betweenCentralities[i];
				largestNodeID = i;
			}
		}
		if (largestIncreasedBC < 0.5) break;// no nodes are found
		isFoundforGBC[largestNodeID] = true;
		foundNodes.push_back(largestNodeID);
		totalGBC += largestIncreasedBC;

		if( (k+1)%10 == 0 ) foundTimes.push_back((clock() - st) / 1000.0);
		//printf("%d th node: %d, covered paths: %.3lf\n", k + 1, 	largestNodeID, largestIncreasedBC);
	}// end for (int k = 0; k < K; ++k)
	for (i = 0; i < n; ++i) // may be less than K nodes are chosen
	{
		if (foundNodes.size() >= K) break;
		if (isFoundforGBC[i] == true) continue;
		foundNodes.push_back(i); // arbitraries 
	}

	totalGBC = (totalGBC / L);
	for (i = 0; i < L; ++i) delete []sampleShortestPathsArray[i];
	return totalGBC; // normalized GBC
}
/**/

double Graph::rankByGroupBCWithSamplePairs(vector<PAIR>& randomPairs, vector<int>& foundNodes, int K, vector<double>& foundTimes)
{// time complexity: O(L m + KLn)
	assert(K > 0);
	foundNodes.clear();
	foundNodes.reserve(K);
	foundTimes.clear();
	int L = randomPairs.size();
	clock_t st = clock();

	int i, j, k, l;
	bool debug = true;
	nonUniformArray<int>  successorsAdjLists; // record successors in the shortest paths
	vector<int> degrees; // degree of each node 
	degrees.reserve(n);
	for (i = 0; i < n; ++i) degrees.push_back(adjLists.size(i));
	successorsAdjLists.allocateMemory(degrees);

	assert(L <= maxSampledSources && n <= MAX_ROWS_OneCC);
	//nonUniformArray<int>* sourcesSuccessorsAdjList[maxSampledSources];
	//int  *sourcesNumShortestPathsBefore[maxSampledSources];
	for (i = 0; i < L; ++i)
	{
		sourcesSuccessorsAdjList[i] = new nonUniformArraySmaller<int>;
		sourcesNumShortestPathsBefore[i] = new int[n];
	}
	for (i = 0; i < L; ++i)
		for (j = 0; j < n; ++j)
			sourcesNumShortestPathsBefore[i][j] = 0;
	bool isOnShortestPathsToDest[MAX_ROWS];
	bool isOnPath;

	int src;
	int u, v, w;
	// the number of nodes that their shortest distances to src have been found
	int sz;
	int* adj;

	double totalGBC = 0;
	double largestIncreasedBC;
	int largestNodeID = 0;

	vector<int> nodesInIncreasingOrder;
	nodesInIncreasingOrder.reserve(n);

	//printf("num nodes: %d, numedges: %d\n", n, numOfEdges() *2 );

	int newDistance;
	int dest;
	bool isReachable;
	for (i = 0; i < L; ++i) // find the shortest path DAG from each source
	{
		//if( i%100 == 0) printf("find %d th shortsest path\n", i + 1);
		src = randomPairs[i].s;
		dest = randomPairs[i].t;

		for (j = 0; j < n; ++j) successorsAdjLists.clear(j);
		// find the shortest paths from each node 
		for (j = 0; j < n; ++j) visited[j] = false;
		visited[src] = true;
		distances[src] = 0;
		queue_clear(); // clear the queue
		queue_push(src);
		nodesInIncreasingOrder.clear();
		nodesInIncreasingOrder.push_back(src);
		isReachable = false;
		while (queue_size() > 0) // bfs
		{
			u = queue_pop();
			if (u == dest) break;

			sz = adjLists.size(u);
			adj = &adjLists.access(u, 0);
			newDistance = distances[u] + 1;
			for (j = 0; j < sz; ++j)
			{
				v = adj[j];
				if (false == visited[v])
				{
					visited[v] = true;
					distances[v] = newDistance;
					queue_push(v);
					nodesInIncreasingOrder.push_back(v);
					if (v == dest)isReachable = true;
				}
				if (distances[v] == newDistance)
					successorsAdjLists.push_back(u, v);
			}
		}

		degrees.clear();
		for (j = 0; j < n; ++j) degrees.push_back(successorsAdjLists.size(j));
		sourcesSuccessorsAdjList[i]->allocateMemory(degrees);

		if (false == isReachable) continue; 

		for (j = 0; j < n; ++j) isOnShortestPathsToDest[j] = false;
		isOnShortestPathsToDest[dest] = true;

		for (j = nodesInIncreasingOrder.size() - 1; j >= 0; --j)
		{
			isOnPath = false;
			v = nodesInIncreasingOrder[j];
			sz = successorsAdjLists.size(v);
			adj = &successorsAdjLists.access(v, 0);
			for (k = 0; k < sz; ++k) // 
			{
				w = adj[k];
				if (true == isOnShortestPathsToDest[w])
				{ // some child is on a shortest path to t
					isOnPath = true;
					sourcesSuccessorsAdjList[i]->push_back(v, w);
				}
			}
			if (true == isOnPath) isOnShortestPathsToDest[v] = true;	
		}

		for (j = 0; j < n; ++j) visited[j] = false;
		visited[src] = true;
		queue_clear(); // clear the queue
		queue_push(src);
		sourcesNumShortestPathsBefore[i][src] = 1;
		while (queue_size() > 0) // bfs
		{
			u = queue_pop();
			sz = sourcesSuccessorsAdjList[i]->size(u);
			adj = &sourcesSuccessorsAdjList[i]->access(u, 0);
			for (j = 0; j < sz; ++j)
			{
				v = adj[j];
				if (false == visited[v])
				{
					visited[v] = true;
					queue_push(v);
				}
				sourcesNumShortestPathsBefore[i][v] += sourcesNumShortestPathsBefore[i][u];
			}
		}
	}
	foundTimes.push_back( (clock()-st)/1000.0);

	// find top-k nodes with the largest GBC from the L pairs
	for (i = 0; i < n; ++i) isFoundforGBC[i] = false;
	for (k = 0; k < K; ++k)
	{
		for (i = 0; i < n; ++i) betweenCentralities[i] = 0;
		for (i = 0; i < L; ++i) // for each pair of nodes
		{
			src = randomPairs[i].s;
			if (true == isFoundforGBC[src])
				continue; // this source node is already in the set of found nodes
			
			dest = randomPairs[i].t;
			if (sourcesNumShortestPathsBefore[i][dest] == 0) continue;

			for (j = 0; j < n; ++j) numShortestPathsAfter[j] = 0;
			numShortestPathsAfter[src] = 1;
			queue_clear(); // clear the queue
			queue_push(src);
			nodesInIncreasingOrder.clear();
			nodesInIncreasingOrder.push_back(src);
			for (j = 0; j < n; ++j) visited[j] = false;
			visited[src] = true;
			for (j = 0; j < foundNodes.size(); ++j) visited[foundNodes[j]] = true;
			while (queue_size() > 0) // bfs
			{
				u = queue_pop();
				sz = sourcesSuccessorsAdjList[i]->size(u);
				adj = &sourcesSuccessorsAdjList[i]->access(u, 0);

				for (j = 0; j < sz; ++j)
				{
					v = adj[j];
					if (false == visited[v])
					{
						visited[v] = true;
						queue_push(v);
						numShortestPathsAfter[v] += numShortestPathsAfter[u];
						nodesInIncreasingOrder.push_back(v);
					}
					else if (false == isFoundforGBC[v])
						numShortestPathsAfter[v] += numShortestPathsAfter[u];
				}
			}
			if( numShortestPathsAfter[dest] == 0 ) continue;// not reachable already

			
			for (j = 0; j < n; ++j) isOnShortestPathsToDest[j] = false;
			isOnShortestPathsToDest[dest] = true;
			for (j = 0; j < n; ++j) currrentBetweenCentralities[j] = 0;
			for (j = nodesInIncreasingOrder.size() - 1; j >= 0; --j)
			{
				isOnPath = false;
				v = nodesInIncreasingOrder[j];
				sz = sourcesSuccessorsAdjList[i]->size(v);
				adj = &sourcesSuccessorsAdjList[i]->access(v, 0);
				for (l = 0; l < sz; ++l)
				{
					w = adj[l];
					if (true == isFoundforGBC[w]) continue;
					if (true == isOnShortestPathsToDest[w]) {
						currrentBetweenCentralities[v] += numShortestPathsAfter[v]
							* currrentBetweenCentralities[w] / numShortestPathsAfter[w];
						isOnPath = true;
					}	
				}
				if( v == dest) currrentBetweenCentralities[v] += 
					     (double)numShortestPathsAfter[v] / sourcesNumShortestPathsBefore[i][v];
				if (true == isOnPath) isOnShortestPathsToDest[v] = true;
			}
			//currrentBetweenCentralities[src]--;
			
			for (j = 0; j < n; ++j)
				betweenCentralities[j] += currrentBetweenCentralities[j];
		}// end for(i = 0; i < L; ++i)
		largestIncreasedBC = 0;
		for (j = 0; j < n; ++j)
		{
			//if(betweenCentralities[j] > 0) assert(isFoundforGBC[j] == false);
			if (betweenCentralities[j] > largestIncreasedBC)
			{
				largestIncreasedBC = betweenCentralities[j];
				largestNodeID = j;
			}
		}

		if (largestIncreasedBC <= 0.0001) break; // no nodes are found

		isFoundforGBC[largestNodeID] = true;
		foundNodes.push_back(largestNodeID);
		totalGBC += largestIncreasedBC;
		if( (k+1)%10 == 0 ) foundTimes.push_back((clock() - st) / 1000.0);
	    //printf("%d th node: %d, increased BC: %.3lf\n", k + 1, 	largestNodeID, largestIncreasedBC);
	}// end for (int k = 0; k < K; ++k)
	for (i = 0; i < n; ++i) // may be less than K nodes are chosen
	{
		if (foundNodes.size() >= K) break;
		if (isFoundforGBC[i] == true) continue;
		foundNodes.push_back(i); // arbitraries 
	}

	totalGBC = totalGBC / L;
	for (i = 0; i < L; ++i)
	{
		delete sourcesSuccessorsAdjList[i];
		delete []sourcesNumShortestPathsBefore[i];
		sourcesSuccessorsAdjList[i] = NULL;
		sourcesNumShortestPathsBefore[i] = NULL;
	}

	return totalGBC;
}/**/


void Graph::recursiveCalGBC(int v)
{
	int w;
	int sz;
	int* adj;
	
	sz = shortestPathsAdjLists->size(v);
	adj = &shortestPathsAdjLists->access(v, 0);

	visited[v] = true;
	for (int i = 0; i< sz; ++i)
	{
		w = adj[i];
		if (true == isFoundforGBC[ w ]) continue;
		if (false == visited[w])
			recursiveCalGBC(w);
		currrentBetweenCentralities[v] += oneSourceNumShortestPathsAfter[v]
			* currrentBetweenCentralities[w] / oneSourceNumShortestPathsAfter[w];
	}
	currrentBetweenCentralities[v] += (double)oneSourceNumShortestPathsAfter[v] / oneSourceNumShortestPathsBefore[v];
}


/*void Graph::findReachNodeByPostOrderDFS(int u, nonUniformArray<int>& adjList, vector<int>& nodesInOrder)
{
	int sz = adjList.size(u);
	int* adj = &adjList.access(u, 0);
	int v;
	for (int j = 0; j < sz; ++j)
	{
		v = adj[j];
		if (false == visited[v]) {
			visited[v] = true;
			findReachNodeByPostOrderDFS(v, adjList, nodesInOrder);
		}
	}
	nodesInOrder.push_back(u);
}

double Graph::rankByGroupBCWithSampleSourcesMemoryFriendlyDFSPostOrder(vector<int>& sourceNodes, vector<int>& foundNodes, int K, vector<double>& foundTimes)
{
	// time complexity: L m + KL 1.5n
	assert(K > 0);
	foundNodes.clear();
	foundNodes.reserve(K);
	foundTimes.clear();
	int L = sourceNodes.size();

	int i, j, k, l;
	bool debug = false;
	nonUniformArray<int>  successorsAdjLists; // record successors in the shortest paths
	vector<int> degrees; // degree of each node 
	degrees.reserve(n);
	for (i = 0; i < n; ++i) degrees.push_back(adjLists.size(i));
	successorsAdjLists.allocateMemory(degrees);

	assert(L <= maxSampledSources && n <= MAX_ROWS_OneCC);
	int* sourcesOut2InIndex = new int[n];
	for (i = 0; i < L; ++i) {
		sourcesSuccessorsAdjList[i] = new nonUniformArraySmaller<int>;
		sourcesNumShortestPathsBefore[i] = new int[n];
		sourcesNumShortestPathsAfter[i] = new int[n];
		sourcesIn2OutIndex[i] = new int[n];
		sourcesPredecessorsAdjList[i] = new nonUniformArraySmaller<int>;
		sourcesGBC[i] = new double[n];
	}
	for (i = 0; i < L; ++i)
		for (j = 0; j < n; ++j)
			sourcesIn2OutIndex[i][j] = -1;

	vector<int> nodesInIncreasingOrder;
	nodesInIncreasingOrder.reserve(n);

	int decreaseDelta[MAX_ROWS];
	for (i = 0; i < L; ++i)
		for (j = 0; j < n; ++j)
			sourcesNumShortestPathsBefore[i][j] = 0;

	int numDecendantsOfOneNode;

	int src;
	int u, v, w;
	// the number of nodes that their shortest distances to src have been found
	int sz;
	int* adj;

	double totalGBC = 0;
	double largestIncreasedBC;
	int outLargestNodeID = 0;

	//printf("num nodes: %d, numedges: %d\n", n, numOfEdges() *2 );
	clock_t st = clock();

	double numEdgesInDAG = 0;
	int newDistance;
	int numReachNodesEachSource;
	double sumGBC;

	double maxMemoryDis = 0;
	const int maxDis = 320;
	const int scale = 32;
	const int numCount = maxDis / scale;
	double disCount[numCount];

	int curDis;
	for (i = 0; i < numCount; ++i) disCount[i] = 0;
	double totalNonZeros = 0;

	int* sourceNodeIDs = new int[L];
	for (i = 0; i < L; ++i) // find the shortest path DAG from each source
	{
		if (i % 10 == 0) printf("finding the %dth shortest-path DAG...\n", i + 1);
		src = sourceNodes[i];
		oneSourceNumShortestPathsBefore = sourcesNumShortestPathsBefore[i];
		oneSourceNumShortestPathsAfter = sourcesNumShortestPathsAfter[i];

		oneSourceNumShortestPathsBefore[src] = 1;
		for (j = 0; j < n; ++j) successorsAdjLists.clear(j);
		// find the shortest paths from each node 
		for (j = 0; j < n; ++j) visited[j] = false;
		visited[src] = true;
		distances[src] = 0;
		queue_clear(); // clear the queue
		queue_push(src);
		while (queue_size() > 0) // bfs, find the shorest path from src to other nodes
		{
			u = queue_pop();
			sz = adjLists.size(u);
			adj = adjLists.startAddress(u);
			newDistance = distances[u] + 1;
			for (j = 0; j < sz; ++j)
			{
				v = adj[j];
				if (false == visited[v])
				{
					visited[v] = true;
					distances[v] = newDistance;
					queue_push(v);
				}
				if (distances[v] == newDistance) {
					successorsAdjLists.push_back(u, v); // create successorsAdjLists in the shortest-path DAG
					oneSourceNumShortestPathsBefore[v] += oneSourceNumShortestPathsBefore[u];
				}
			}
		}

		nodesInIncreasingOrder.clear();
		for (j = 0; j < n; ++j) visited[j] = false;
		visited[src] = true;
		findReachNodeByPostOrderDFS(src, successorsAdjLists, nodesInIncreasingOrder);

		numReachNodesEachSource = nodesInIncreasingOrder.size(); // reorder the nodes in the shortest-path DAG
		assert(numReachNodesEachSource <= n);
		sourceNodeIDs[i] = numReachNodesEachSource-1;
		for (j = 0; j < numReachNodesEachSource; ++j) {
			v = nodesInIncreasingOrder[j]; // reorder the nodes by the post-order DFS order
			sourcesIn2OutIndex[i][j] = v;
			sourcesOut2InIndex[v] = j;
		}

		degrees.clear();
		for (j = 0; j < numReachNodesEachSource; ++j)
			degrees.push_back(successorsAdjLists.size(nodesInIncreasingOrder[j]));
		shortestPathsAdjLists = sourcesSuccessorsAdjList[i];
		shortestPathsAdjLists->allocateMemory(degrees);
		degrees.clear();
		degrees.assign(n, 0); // cal degrees for predecessor lists
		for (j = 0; j < numReachNodesEachSource; ++j)
		{
			v = nodesInIncreasingOrder[j];
			sz = successorsAdjLists.size(v);
			adj = successorsAdjLists.startAddress(v);
			for (k = 0; k < sz; ++k) {
				w = sourcesOut2InIndex[adj[k]];
				shortestPathsAdjLists->push_back(j, w);
				degrees[w]++; // count the degrees for the predecessor lists
			}
			oneSourceNumShortestPathsAfter[j] = oneSourceNumShortestPathsBefore[v];
			numEdgesInDAG += shortestPathsAdjLists->size(j);
		}
		for (j = 0; j < numReachNodesEachSource; ++j) // make them same
			oneSourceNumShortestPathsBefore[j] = oneSourceNumShortestPathsAfter[j];

		shortestPathsPredecessorAdjLists = sourcesPredecessorsAdjList[i];
		shortestPathsPredecessorAdjLists->allocateMemory(degrees); // create predecessor lists
		for (j = 0; j < numReachNodesEachSource; ++j)
		{
			sz = shortestPathsAdjLists->size(j);
			adj = shortestPathsAdjLists->startAddress(j);
			for (k = 0; k < sz; ++k) 
				shortestPathsPredecessorAdjLists->push_back(adj[k], j); // node j is a predecessor of adj[k]
		}

	}
	successorsAdjLists.freeMemory();
	delete[]sourcesOut2InIndex;
	sourcesOut2InIndex = NULL;

	numEdgesInDAG /= L;
	clock_t ft1 = clock();
	printf("n: %d, numEdgesInDAG:%.0lf, ratio: %.2lf, time for finding shortest paths: %.3lf,\n", n, numEdgesInDAG, numEdgesInDAG / n, (ft1 - st) / 1000.0);
	foundTimes.push_back((ft1 - st) / 1000.0);

	for (i = 0; i < L; ++i) {
		oneSourceGBC = sourcesGBC[i];// calculate the betweenness centrality of each node
		for (j = 0; j < n; ++j) oneSourceGBC[j] = 0;
		oneSourceNumShortestPathsBefore = sourcesNumShortestPathsBefore[i];
		shortestPathsAdjLists = sourcesSuccessorsAdjList[i];
		for (j = 0; j < n; ++j)// notice that the nodes have been sorted in the post-order DFS
		{
			if (0 == oneSourceNumShortestPathsBefore[j]) continue;

			sumGBC = 0;
			sz = shortestPathsAdjLists->size(j);
			adj = shortestPathsAdjLists->startAddress(j);
			for (l = 0; l < sz; ++l)
			{
				w = adj[l];
				sumGBC += oneSourceNumShortestPathsBefore[j] * oneSourceGBC[w] / oneSourceNumShortestPathsBefore[w];
			}
			oneSourceGBC[j] = sumGBC + 1;
		}
		oneSourceGBC[ sourceNodeIDs[i] ]--; // the last node is the source node
	}

	int innerLargestID;
	int* In2OutIndex;
	// find top-k nodes with the largest GBC from the L sources
	for (i = 0; i < n; ++i) isFoundforGBC[i] = false;

	vector<int> nodesNeedUpdateGBC;
	nodesNeedUpdateGBC.reserve(n);
	int* allSourcesNumReachNodes = new int[L];
	for (i = 0; i < L; ++i)allSourcesNumReachNodes[i] = n;

	int halfSize, numNodesNeedUpdate;
	int tmp, reverseIndex;
	int sourceID;
	for (k = 0; k < K; ++k)
	{
		for (i = 0; i < n; ++i) betweenCentralities[i] = 0;  // outer indices
		for (i = 0; i < L; ++i) // calculate the increased GBC of each node
		{
			if (true == isFoundforGBC[sourceNodes[i]])
				continue; // this source node is already in the set of found nodes
			In2OutIndex = sourcesIn2OutIndex[i];
			oneSourceGBC = sourcesGBC[i];
			for (j = 0; j < n; ++j)
				betweenCentralities[In2OutIndex[j]] += oneSourceGBC[j];
		}

		largestIncreasedBC = 0;
		outLargestNodeID = -1;
		for (j = 0; j < n; ++j)
			if (betweenCentralities[j] > largestIncreasedBC)
			{
				largestIncreasedBC = betweenCentralities[j];
				outLargestNodeID = j;
			}
		if (-1 == outLargestNodeID) break; // no more nodes 

		//if ((k + 1) % 10 == 0)
		printf("find the %d th node: %d, increased GBC: %.0lf, used total time: %.2lf, delta Time: %.3lf \n", k + 1, outLargestNodeID, largestIncreasedBC, (clock() - st) / 1000.0, (clock() - ft1) / 1000.0);
		ft1 = clock();
		if ((k + 1) % 10 == 0) foundTimes.push_back((ft1 - st) / 1000.0);
		foundNodes.push_back(outLargestNodeID);
		totalGBC += largestIncreasedBC;

		if (k == K - 1) {
			isFoundforGBC[outLargestNodeID] = true;
			break; // is the last node, no need to update GBC
		}

		for (j = 0; j < n; ++j) decreaseDelta[j] = 0;
		for (j = 0; j < n; ++j) visited[j] = false;
		for (i = 0; i < L; ++i) // updating
		{
			if (true == isFoundforGBC[sourceNodes[i]])
				continue; // this source node is already in the set of found nodes
			// update the number of shortest paths of each node
			shortestPathsAdjLists = sourcesSuccessorsAdjList[i];
			shortestPathsPredecessorAdjLists = sourcesPredecessorsAdjList[i];
			oneSourceNumShortestPathsBefore = sourcesNumShortestPathsBefore[i];
			oneSourceNumShortestPathsAfter = sourcesNumShortestPathsAfter[i];
			oneSourceGBC = sourcesGBC[i];

			In2OutIndex = sourcesIn2OutIndex[i];
			for (j = 0; j < n; ++j) // find the inner node ID of the found node in this shorest-path DAG
				if (In2OutIndex[j] == outLargestNodeID) break;
			assert(j != n);
			innerLargestID = j;

			decreaseDelta[innerLargestID] = oneSourceNumShortestPathsAfter[innerLargestID];
			if (0 == decreaseDelta[innerLargestID]) continue; // no need to update

			if (sourceNodeIDs[i] == innerLargestID) { // is the source node
				for (j = 0; j < n; ++j) {
					oneSourceNumShortestPathsAfter[j] = 0;
					oneSourceGBC[j] = 0;
				}
				continue;
			}
			visited[innerLargestID] = true;
			nodesNeedUpdateGBC.clear();
			nodesNeedUpdateGBC.push_back(innerLargestID);// record the set of reachable nodes in this shorest-path DAG			

			queue_clear(); // clear the queue
			queue_push(innerLargestID); // the kth found node
			while (queue_size() > 0) // bfs, decrease the number of shorestest paths
			{ // find the set of reachable nodes in this shorest-path DAG
				u = queue_pop();
				sz = shortestPathsAdjLists->size(u);
				adj = &shortestPathsAdjLists->access(u, 0);
				for (j = 0; j < sz; ++j)
				{
					v = adj[j];
					if (0 == oneSourceNumShortestPathsAfter[v]) continue;
					if (false == visited[v]) {
						queue_push(v);
						visited[v] = true;
						nodesNeedUpdateGBC.push_back(v);
					}
					decreaseDelta[v] += decreaseDelta[u]; // the number of decreased shortest paths
				}
			}
			numDecendantsOfOneNode = nodesNeedUpdateGBC.size();

			if (numDecendantsOfOneNode >= allSourcesNumReachNodes[i] / 10) {
				//if( k<0){
				for (j = 0; j < numDecendantsOfOneNode; ++j) { // decrease the number of shoretest paths
					v = nodesNeedUpdateGBC[j];
					oneSourceNumShortestPathsAfter[v] -= decreaseDelta[v];
					decreaseDelta[v] = 0;
					visited[v] = false;
					if (0 == oneSourceNumShortestPathsAfter[v])  allSourcesNumReachNodes[i]--;
				}
				// calculate the betweenness centrality of each node
				for (j = 0; j < n; ++j) oneSourceGBC[j] = 0;
				for (j = 0; j < n; ++j)// notice that the nodes have been sorted in the DFS post-order
				{
					if (0 == oneSourceNumShortestPathsAfter[j]) continue;
					sumGBC = 0;
					sz = shortestPathsAdjLists->size(j);
					adj = shortestPathsAdjLists->startAddress(j);
					for (l = 0; l < sz; ++l)
					{
						w = adj[l];
						if (0 == oneSourceNumShortestPathsAfter[w]) continue;
						sumGBC += oneSourceNumShortestPathsAfter[j] * oneSourceGBC[w] / oneSourceNumShortestPathsAfter[w];
					}
					oneSourceGBC[j] = sumGBC + double(oneSourceNumShortestPathsAfter[j]) / oneSourceNumShortestPathsBefore[j];
				}
				if (innerLargestID != sourceNodeIDs[i]) oneSourceGBC[sourceNodeIDs[i]]--; //   source node

				continue;//
			}

			halfSize = numDecendantsOfOneNode / 2; // revserse the order. nodes are in increasing order
			for (j = 0, reverseIndex = numDecendantsOfOneNode - 1; j < halfSize; ++j, --reverseIndex) {
				tmp = nodesNeedUpdateGBC[j];
				nodesNeedUpdateGBC[j] = nodesNeedUpdateGBC[reverseIndex];
				nodesNeedUpdateGBC[reverseIndex] = tmp;
			}

			queue_clear(); // find the rest nodes that need update GBC
			for (j = 0; j < numDecendantsOfOneNode; ++j)
				queue_push(nodesNeedUpdateGBC[j]);
			while (queue_size() > 0) // backward BFS
			{
				v = queue_pop();
				sz = shortestPathsPredecessorAdjLists->size(v);
				adj = shortestPathsPredecessorAdjLists->startAddress(v);
				for (j = 0; j < sz; ++j)
				{
					u = adj[j]; // each predecessor
					if (0 == oneSourceNumShortestPathsAfter[u]) continue;// this predecessor is not reachable from src
					if (false == visited[u]) {
						queue_push(u);
						visited[u] = true;
						nodesNeedUpdateGBC.push_back(u);
					}
				}
			}

			for (j = 0; j < numDecendantsOfOneNode; ++j) { // decrease the number of shoretest paths
				v = nodesNeedUpdateGBC[j];
				oneSourceNumShortestPathsAfter[v] -= decreaseDelta[v];
				decreaseDelta[v] = 0;
				if (0 == oneSourceNumShortestPathsAfter[v])  allSourcesNumReachNodes[i]--;
			}

			numNodesNeedUpdate = nodesNeedUpdateGBC.size();
			for (j = 0; j < numNodesNeedUpdate; ++j) // nodes are partially sorted in the DFS post-order
			{
				v = nodesNeedUpdateGBC[j];
				if (0 == oneSourceNumShortestPathsAfter[v]) {
					oneSourceGBC[v] = 0;
					continue;
				}
				sumGBC = 0;
				sz = shortestPathsAdjLists->size(v);
				adj = shortestPathsAdjLists->startAddress(v);
				for (l = 0; l < sz; ++l)
				{
					w = adj[l];
					if (0 == oneSourceNumShortestPathsAfter[w]) continue;
					sumGBC += oneSourceNumShortestPathsAfter[v] * oneSourceGBC[w] / oneSourceNumShortestPathsAfter[w];
				}
				oneSourceGBC[v] = sumGBC + double(oneSourceNumShortestPathsAfter[v]) / oneSourceNumShortestPathsBefore[v];
			}
			if (innerLargestID != sourceNodeIDs[i])
				oneSourceGBC[ sourceNodeIDs[i] ]--; //  the source nodes must be updated

			for (j = 0; j < nodesNeedUpdateGBC.size(); ++j)
				visited[nodesNeedUpdateGBC[j]] = false;
		} //for (i = 0; i < L; ++i) // updating

		isFoundforGBC[outLargestNodeID] = true;
	}// end for (int k = 0; k < K; ++k)
	for (i = 0; i < n; ++i) // there may be less than K nodes in the set of found nodes
	{
		if (foundNodes.size() >= K) break;
		if (false == isFoundforGBC[i]) foundNodes.push_back(i);
	}

	double avgDis = 0;
	totalNonZeros = 1;
	for (i = 0; i < numCount; ++i) {
		disCount[i] /= totalNonZeros;
		avgDis += disCount[i] * (i * scale);
	}

	for (i = 0; i < numCount; ++i)
		printf("interval: [%d, %d): %.6lf %%\n", i * scale, (i + 1) * scale, 100 * disCount[i]);
	printf("\navgDis: %.6lf\n\n", avgDis);

	totalGBC = (totalGBC / L);
	totalGBC /= double(n - 1.0); // normalized GBC

	for (i = 0; i < L; ++i) {
		delete sourcesSuccessorsAdjList[i];
		sourcesSuccessorsAdjList[i] = NULL;
		delete[] sourcesNumShortestPathsBefore[i];
		sourcesNumShortestPathsBefore[i] = NULL;
		delete[]sourcesNumShortestPathsAfter[i];
		sourcesNumShortestPathsAfter[i] = NULL;
		delete[]sourcesIn2OutIndex[i];
		sourcesIn2OutIndex[i] = NULL;
		delete sourcesPredecessorsAdjList[i];
		sourcesPredecessorsAdjList[i] = NULL;
		delete[]sourcesGBC[i];
		sourcesGBC[i] = NULL;
	}
	delete[]allSourcesNumReachNodes;
	delete[]sourceNodeIDs;

	return totalGBC;
}
/**/

double Graph::rankByGroupBCWithSampleSourcesMemoryFriendlyBFS(vector<int>& sourceNodes, vector<int>& foundNodes, int K, vector<double>& foundTimes)
{
	// time complexity: L m + KL 1.5n
	assert(K > 0);
	foundNodes.clear();
	foundNodes.reserve(K);
	foundTimes.clear();
	int L = sourceNodes.size();

	int i, j, k, l;
	bool debug = false;
	nonUniformArray<int>  successorsAdjLists; // record successors in the shortest paths
	vector<int> degrees; // degree of each node 
	degrees.reserve(n);
	for (i = 0; i < n; ++i) degrees.push_back(adjLists.size(i));
	successorsAdjLists.allocateMemory(degrees);

	assert(L <= maxSampledSources && n <= MAX_ROWS_OneCC);
	int* sourcesOut2InIndex = new int[n];
	for (i = 0; i < L; ++i) {
		sourcesSuccessorsAdjList[i] = new nonUniformArraySmaller<int>;
		sourcesNumShortestPathsBefore[i] = new int[n];
		sourcesNumShortestPathsAfter[i] = new int[n];
		sourcesIn2OutIndex[i] = new int[n];
		sourcesPredecessorsAdjList[i] = new nonUniformArraySmaller<int>;
		sourcesGBC[i] = new double[n];
	}
	for (i = 0; i < L; ++i)
		for (j = 0; j < n; ++j)
			sourcesIn2OutIndex[i][j] = -1;

	vector<int> nodesInIncreasingOrder;
	nodesInIncreasingOrder.reserve(n);

	int decreaseDelta[MAX_ROWS];
	for (i = 0; i < L; ++i)
		for (j = 0; j < n; ++j)
			sourcesNumShortestPathsBefore[i][j] = 0;
	
	int numDecendantsOfOneNode;
	int src;
	int u, v, w;
	// the number of nodes that their shortest distances to src have been found
	int sz;
	int* adj;

	double totalGBC = 0;
	double largestIncreasedBC;
	int outLargestNodeID = 0;

	//printf("num nodes: %d, numedges: %d\n", n, numOfEdges() *2 );
	clock_t st = clock();

	double numEdgesInDAG = 0;
	int newDistance;
	int numReachNodesEachSource;
	double sumGBC;

	double maxMemoryDis = 0;
	const int maxDis = 320;
	const int scale = 32;
	const int numCount = maxDis / scale;
	double disCount[numCount];
	int curDis;
	for (i = 0; i < numCount; ++i) disCount[i] = 0;
	double totalNonZeros = 0;

	int* numNodesDifferentSources = new int[L];
	int* allSourcesRestReachNodes = new int[L];

	for (i = 0; i < L; ++i) // find the shortest path DAG from each source
	{
		//if (i % 10 == 0) printf("finding the %dth shortest-path DAG...\n", i + 1);
		src = sourceNodes[i];

		oneSourceNumShortestPathsBefore = sourcesNumShortestPathsBefore[i];
		oneSourceNumShortestPathsAfter = sourcesNumShortestPathsAfter[i];

		oneSourceNumShortestPathsBefore[src] = 1;
		for (j = 0; j < n; ++j) successorsAdjLists.clear(j);
		// find the shortest paths from each node 
		for (j = 0; j < n; ++j) visited[j] = false;
		visited[src] = true;
		distances[src] = 0;
		queue_clear(); // clear the queue
		queue_push(src);
		nodesInIncreasingOrder.clear();
		nodesInIncreasingOrder.push_back(src);
		while (queue_size() > 0) // bfs, find the shorest path from src to other nodes
		{
			u = queue_pop();
			sz = adjLists.size(u);
			adj = adjLists.startAddress(u);
			newDistance = distances[u] + 1;
			for (j = 0; j < sz; ++j)
			{
				v = adj[j];
				if (false == visited[v])
				{
					visited[v] = true;
					distances[v] = newDistance;
					queue_push(v);
					nodesInIncreasingOrder.push_back(v);
				}
				if (distances[v] == newDistance) {
					successorsAdjLists.push_back(u, v); // create successorsAdjLists in the shortest-path DAG
					oneSourceNumShortestPathsBefore[v] += oneSourceNumShortestPathsBefore[u];
				}
			}
		}

		numReachNodesEachSource = nodesInIncreasingOrder.size(); // reorder the nodes in the shortest-path DAG
		numNodesDifferentSources[i] = numReachNodesEachSource;
		allSourcesRestReachNodes[i] = numReachNodesEachSource;
		for (j = 0; j < numReachNodesEachSource; ++j) {
			v = nodesInIncreasingOrder[j]; // reorder the nodes by the BFS order
			sourcesIn2OutIndex[i][j] = v;
			sourcesOut2InIndex[v] = j;
		}

		degrees.clear();
		for (j = 0; j < numReachNodesEachSource; ++j)
			degrees.push_back(successorsAdjLists.size(nodesInIncreasingOrder[j]));
		shortestPathsAdjLists = sourcesSuccessorsAdjList[i];
		shortestPathsAdjLists->allocateMemory(degrees);
		degrees.clear();
		degrees.assign(n, 0); // cal degrees for predecessor lists
		for (j = 0; j < numReachNodesEachSource; ++j)
		{
			v = nodesInIncreasingOrder[j];
			sz = successorsAdjLists.size(v);
			adj = successorsAdjLists.startAddress(v);
			for (k = 0; k < sz; ++k) {
				w = sourcesOut2InIndex[adj[k]];
				shortestPathsAdjLists->push_back(j, w);
				degrees[w]++; // count the degrees for the predecessor lists
			}
			oneSourceNumShortestPathsAfter[j] = oneSourceNumShortestPathsBefore[v];
			numEdgesInDAG += shortestPathsAdjLists->size(j);
		}
		for (j = 0; j < numReachNodesEachSource; ++j) // make them same
			oneSourceNumShortestPathsBefore[j] = oneSourceNumShortestPathsAfter[j];
		
		shortestPathsPredecessorAdjLists = sourcesPredecessorsAdjList[i];
		shortestPathsPredecessorAdjLists->allocateMemory(degrees); // create predecessor lists
		for (j = 0; j < numReachNodesEachSource; ++j)
		{	
			sz = shortestPathsAdjLists->size(j);
			adj = shortestPathsAdjLists->startAddress(j);
			for (k = 0; k < sz; ++k)
				shortestPathsPredecessorAdjLists->push_back(adj[k], j); // node j is a predecessor of adj[k]
		}
		/*printf("src: %d\n", src);
		printf("nodesInIncreasingOrder:");
		for (j = 0; j < nodesInIncreasingOrder.size(); ++j) printf("%d, ", nodesInIncreasingOrder[j]);
		printf("\n");
		shortestPathsAdjLists->printArray();
		printf("\n");
		shortestPathsPredecessorAdjLists->printArray();/**/
	}
	successorsAdjLists.freeMemory();
	delete[]sourcesOut2InIndex;
	sourcesOut2InIndex = NULL;

	numEdgesInDAG /= L;
	clock_t ft1 = clock();
	//printf("n: %d, numEdgesInDAG:%.0lf, ratio: %.2lf, time for finding shortest paths: %.3lf,\n", n, numEdgesInDAG, numEdgesInDAG / n, (ft1 - st) / 1000.0);
	foundTimes.push_back((ft1 - st) / 1000.0);

	for (i = 0; i < L; ++i) {
		oneSourceGBC = sourcesGBC[i];// calculate the betweenness centrality of each node
		for (j = 0; j < n; ++j) oneSourceGBC[j] = 0;
		oneSourceNumShortestPathsBefore = sourcesNumShortestPathsBefore[i];
		shortestPathsAdjLists = sourcesSuccessorsAdjList[i];
		for (j = numNodesDifferentSources[i] - 1; j >= 0; --j)// notice that the nodes have been sorted in the BFS order
		{
			if (0 == oneSourceNumShortestPathsBefore[j]) continue;
			sumGBC = 0;
			sz = shortestPathsAdjLists->size(j);
			adj = shortestPathsAdjLists->startAddress(j);
			for (l = 0; l < sz; ++l)
			{
				w = adj[l];
				sumGBC += oneSourceGBC[w] * oneSourceNumShortestPathsBefore[j] / oneSourceNumShortestPathsBefore[w];
			}
			oneSourceGBC[j] = sumGBC + 1;
		}
		oneSourceGBC[0]--; // node 0 is the source node
	}

	int innerLargestID;
	int* In2OutIndex;
	// find top-k nodes with the largest GBC from the L sources
	for (i = 0; i < n; ++i) isFoundforGBC[i] = false;

	vector<int> nodesNeedUpdateGBC;
	nodesNeedUpdateGBC.reserve(n);

	int halfSize, numNodesNeedUpdate;
	int tmp, reverseIndex;
	for (k = 0; k < K; ++k)
	{
		for (i = 0; i < n; ++i) betweenCentralities[i] = 0;  // outer indices
		for (i = 0; i < L; ++i) // calculate the increased GBC of each node
		{
			if (true == isFoundforGBC[sourceNodes[i]])
				continue; // this source node is already in the set of found nodes
			In2OutIndex = sourcesIn2OutIndex[i];
			oneSourceGBC = sourcesGBC[i];
			numReachNodesEachSource = numNodesDifferentSources[i];
			for (j = 0; j < numReachNodesEachSource; ++j)
				betweenCentralities[ In2OutIndex[j] ] += oneSourceGBC[j];
		}

		largestIncreasedBC = 0;
		outLargestNodeID = -1;
		for (j = 0; j < n; ++j)
			if (betweenCentralities[j] > largestIncreasedBC)
			{
				largestIncreasedBC = betweenCentralities[j];
				outLargestNodeID = j;
			}
		if (-1 == outLargestNodeID) break; // no more nodes 

		foundNodes.push_back(outLargestNodeID);
		totalGBC += largestIncreasedBC;
		if ((k + 1) % 10 == 0) foundTimes.push_back((ft1 - st) / 1000.0);

		//if ((k + 1) % 10 == 0)printf("find the %d th node: %d, increased GBC: %.0lf, used total time: %.2lf, delta Time: %.3lf \n", k + 1, outLargestNodeID, largestIncreasedBC, (clock() - st) / 1000.0, (clock() - ft1) / 1000.0);
		ft1 = clock();

		if (k == K - 1) {
			isFoundforGBC[outLargestNodeID] = true;
			break; // is the last node, no need to update GBC
		}
		
		for (j = 0; j < n; ++j) decreaseDelta[j] = 0;
		for (j = 0; j < n; ++j) visited[j] = false;
		for (i = 0; i < L; ++i) // updating
		{ 
			if (true == isFoundforGBC[sourceNodes[i]])
				continue; // this source node is already in the set of found nodes
			// update the number of shortest paths of each node
			shortestPathsAdjLists = sourcesSuccessorsAdjList[i];
			shortestPathsPredecessorAdjLists = sourcesPredecessorsAdjList[i];
			oneSourceNumShortestPathsBefore = sourcesNumShortestPathsBefore[i];
			oneSourceNumShortestPathsAfter = sourcesNumShortestPathsAfter[i];
			oneSourceGBC = sourcesGBC[i];

			In2OutIndex = sourcesIn2OutIndex[i]; 
			numReachNodesEachSource = numNodesDifferentSources[i];
			for (j = 0; j < numReachNodesEachSource; ++j) // find the inner node ID of the found node in this shorest-path DAG
				if (In2OutIndex[j] == outLargestNodeID) break;
			assert(j != n);
			innerLargestID = j;

			decreaseDelta[innerLargestID] = oneSourceNumShortestPathsAfter[innerLargestID];
			if (0 == decreaseDelta[innerLargestID]) continue; // no need to update

			if (0 == innerLargestID) { // node 0 is the source node
				numReachNodesEachSource =  numNodesDifferentSources[i];
				for (j = 0; j < numReachNodesEachSource; ++j) {
					oneSourceNumShortestPathsAfter[j] = 0;
					oneSourceGBC[j] = 0;
				}
				allSourcesRestReachNodes[i] = 0;
				continue;
			}

			visited[innerLargestID] = true;
			nodesNeedUpdateGBC.clear();
			nodesNeedUpdateGBC.push_back(innerLargestID);// record the set of reachable nodes in this shorest-path DAG			
			queue_clear(); // clear the queue
			queue_push(innerLargestID); // the kth found node
			while (queue_size() > 0) // bfs, decrease the number of shorestest paths
			{ // find the set of reachable nodes in this shorest-path DAG
				u = queue_pop();
				sz = shortestPathsAdjLists->size(u);
				adj = shortestPathsAdjLists->startAddress(u);
				for (j = 0; j < sz; ++j)
				{
					v = adj[j];
					if (0 == oneSourceNumShortestPathsAfter[v]) continue;
					if (false == visited[v]) {
						queue_push(v);
						visited[v] = true;
						nodesNeedUpdateGBC.push_back(v);
					}
					decreaseDelta[v] += decreaseDelta[u]; // the number of decreased shortest paths
				}
			}
			numDecendantsOfOneNode = nodesNeedUpdateGBC.size();

			if (numDecendantsOfOneNode >= allSourcesRestReachNodes[i]/10) {
			//if (numDecendantsOfOneNode >= numNodesDifferentSources[i] / 10) {
			//if( k<0){
				for (j = 0; j < numDecendantsOfOneNode; ++j) { // decrease the number of shoretest paths
					v = nodesNeedUpdateGBC[j];
					oneSourceNumShortestPathsAfter[v] -= decreaseDelta[v];
					decreaseDelta[v] = 0;
					visited[v] = false;
					if (0 == oneSourceNumShortestPathsAfter[v])  allSourcesRestReachNodes[i]--;
				}
				// calculate the betweenness centrality of each node
				numReachNodesEachSource = numNodesDifferentSources[i];
				for (j = 0; j < numReachNodesEachSource; ++j) oneSourceGBC[j] = 0;
				for (j = numReachNodesEachSource - 1; j >= 0; --j)// notice that the nodes have been sorted in the BFS order
				{
					if (0 == oneSourceNumShortestPathsAfter[j]) continue;
					sumGBC = 0;
					sz = shortestPathsAdjLists->size(j);
					adj = shortestPathsAdjLists->startAddress(j);
					for (l = 0; l < sz; ++l)
					{
						w = adj[l];
						if (0 == oneSourceNumShortestPathsAfter[w]) continue;
						sumGBC += oneSourceNumShortestPathsAfter[j] * oneSourceGBC[w] / oneSourceNumShortestPathsAfter[w];
					}
					oneSourceGBC[j] = sumGBC + double(oneSourceNumShortestPathsAfter[j])/ oneSourceNumShortestPathsBefore[j];
				}
				oneSourceGBC[0]--; // node 0 is the source node
				continue;//
			}

			halfSize = numDecendantsOfOneNode / 2; // revserse the order. nodes are in decreasing order
			for (j = 0, reverseIndex = numDecendantsOfOneNode-1; j < halfSize; ++j, --reverseIndex) {
				tmp = nodesNeedUpdateGBC[j];
				nodesNeedUpdateGBC[j] = nodesNeedUpdateGBC[reverseIndex];
				nodesNeedUpdateGBC[reverseIndex] = tmp;
			}

			queue_clear(); // find the rest nodes that need update GBC
			for (j = 0; j < numDecendantsOfOneNode; ++j)
				queue_push(nodesNeedUpdateGBC[j]);
			while (queue_size() > 0) // backward BFS
			{
				v = queue_pop();
				sz = shortestPathsPredecessorAdjLists->size(v);
				adj = shortestPathsPredecessorAdjLists->startAddress(v);
				for (j = 0; j < sz; ++j)
				{
					u = adj[j]; // each predecessor
					if (0 == oneSourceNumShortestPathsAfter[u]) continue;// this predecessor is not reachable from src
					if (false == visited[u]) {
						queue_push(u);
						visited[u] = true;
						nodesNeedUpdateGBC.push_back(u);
					}
				}
			}

			for (j = 0; j < numDecendantsOfOneNode; ++j) { // decrease the number of shoretest paths
				v = nodesNeedUpdateGBC[j];
				oneSourceNumShortestPathsAfter[v] -= decreaseDelta[v];
				decreaseDelta[v] = 0;
				if (0 == oneSourceNumShortestPathsAfter[v])  allSourcesRestReachNodes[i]--;
			}

			numNodesNeedUpdate = nodesNeedUpdateGBC.size();
			for (j = 0; j < numNodesNeedUpdate;++j) // nodes are partially sorted in the BFS order
			{
				v = nodesNeedUpdateGBC[j];
				if (0 == oneSourceNumShortestPathsAfter[v]) {
					oneSourceGBC[v] = 0;
					continue;
				}
				sumGBC = 0;
				sz = shortestPathsAdjLists->size(v);
				adj = shortestPathsAdjLists->startAddress(v);
				for (l = 0; l < sz; ++l)
				{
					w = adj[l];
					if (0 == oneSourceNumShortestPathsAfter[w]) continue;
					sumGBC += oneSourceNumShortestPathsAfter[v] * oneSourceGBC[w] / oneSourceNumShortestPathsAfter[w];
				}
				oneSourceGBC[v] = sumGBC + double(oneSourceNumShortestPathsAfter[v])/ oneSourceNumShortestPathsBefore[v];
			}
			oneSourceGBC[0]--; // node 0 is the source node, the source nodes must be updated

			for (j = 0; j < nodesNeedUpdateGBC.size(); ++j)
				visited[nodesNeedUpdateGBC[j]] = false;
		} //for (i = 0; i < L; ++i) // updating

		isFoundforGBC[outLargestNodeID] = true;
	}// end for (int k = 0; k < K; ++k)
	for (i = 0; i < n; ++i) // there may be less than K nodes in the set of found nodes
	{
		if (foundNodes.size() >= K) break;
		if (false == isFoundforGBC[i]) foundNodes.push_back(i);
	}

	double avgDis = 0;
	totalNonZeros = 1;
	for (i = 0; i < numCount; ++i) {
		disCount[i] /= totalNonZeros;
		avgDis += disCount[i] * (i * scale);
	}

//	for (i = 0; i < numCount; ++i)
		//printf("interval: [%d, %d): %.6lf %%\n", i * scale, (i + 1) * scale, 100 * disCount[i]);printf("\navgDis: %.6lf\n\n", avgDis);

	totalGBC = (totalGBC / L);
	totalGBC /= double(n - 1.0); // normalized GBC

	for (i = 0; i < L; ++i) {
		delete sourcesSuccessorsAdjList[i];
		sourcesSuccessorsAdjList[i] = NULL;
		delete[] sourcesNumShortestPathsBefore[i];
		sourcesNumShortestPathsBefore[i] = NULL;
		delete[]sourcesNumShortestPathsAfter[i];
		sourcesNumShortestPathsAfter[i] = NULL;
		delete[]sourcesIn2OutIndex[i]; 
		sourcesIn2OutIndex[i] = NULL;
		delete sourcesPredecessorsAdjList[i];
		sourcesPredecessorsAdjList[i] = NULL;
		delete[]sourcesGBC[i];
		sourcesGBC[i] = NULL;
	}
	delete[]allSourcesRestReachNodes;
	delete[]numNodesDifferentSources;

	return totalGBC;
}


double Graph::rankByGroupBCWithSampleSourcesMemoryFriendly(vector<int>& sourceNodes, vector<int>& foundNodes, int K, vector<double>& foundTimes)
{
	// time complexity: L m + KL 1.5n
	assert(K > 0);
	foundNodes.clear();
	foundNodes.reserve(K);
	foundTimes.clear();
	int L = sourceNodes.size();

	int i, j, k, l;
	bool debug = false;
	nonUniformArray<int>  successorsAdjLists; // record successors in the shortest paths
	vector<int> degrees; // degree of each node 
	degrees.reserve(n);
	for (i = 0; i < n; ++i) degrees.push_back(adjLists.size(i));
	successorsAdjLists.allocateMemory(degrees);

	assert(L <= maxSampledSources && n <= MAX_ROWS_OneCC);
	int* sourcesOut2InIndex = new int[n];
	for (i = 0; i < L; ++i) {
		sourcesSuccessorsAdjList[i] = new nonUniformArraySmaller<int>;
		sourcesNumShortestPathsBefore[i] = new int[n];
		sourcesNumShortestPathsAfter[i] = new int[n];
		sourcesIn2OutIndex[i] = new int[n];
	}
	for (i = 0; i < L; ++i) 
		for (j = 0; j < n; ++j) 
			sourcesIn2OutIndex[i][j] = -1;

	vector<int> nodesInIncreasingOrder;
	nodesInIncreasingOrder.reserve(n);

	int decreaseDelta[MAX_ROWS];
	for (i = 0; i < L; ++i)
		for (j = 0; j < n; ++j)
			sourcesNumShortestPathsBefore[i][j] = 0;
	int recordNodes[MAX_ROWS + 1000];
	int numRecordNodes;

	int src;
	int u, v, w;
	// the number of nodes that their shortest distances to src have been found
	int sz;
	int* adj;

	double totalGBC = 0;
	double largestIncreasedBC;
	int outLargestNodeID = 0;

	//printf("num nodes: %d, numedges: %d\n", n, numOfEdges() *2 );
	clock_t st = clock();

	double numEdgesInDAG = 0;
	int newDistance;
	int numReachNodesEachSource;
	
	for (i = 0; i < L; ++i) // find the shortest path DAG from each source
	{
		if (i % 10 == 0) printf("finding the %dth shortest-path DAG...\n", i + 1);
		src = sourceNodes[i];

		sourcesNumShortestPathsBefore[i][src] = 1;
		for (j = 0; j < n; ++j) successorsAdjLists.clear(j);
		// find the shortest paths from each node 
		for (j = 0; j < n; ++j) visited[j] = false;
		visited[src] = true;
		distances[src] = 0;
		queue_clear(); // clear the queue
		queue_push(src);
		nodesInIncreasingOrder.clear();
		nodesInIncreasingOrder.push_back(src);
		while (queue_size() > 0) // bfs
		{
			u = queue_pop();

			sz = adjLists.size(u);
			adj = &adjLists.access(u, 0);
			newDistance = distances[u] + 1;
			for (j = 0; j < sz; ++j)
			{
				v = adj[j];
				if (false == visited[v])
				{
					visited[v] = true;
					distances[v] = newDistance;
					queue_push(v);
					nodesInIncreasingOrder.push_back(v);
				}
				if (distances[v] == newDistance) {
					successorsAdjLists.push_back(u, v);
					sourcesNumShortestPathsBefore[i][v] += sourcesNumShortestPathsBefore[i][u];
				}
			}
		}

		numReachNodesEachSource = nodesInIncreasingOrder.size();
		assert(numReachNodesEachSource <= n);
		for (j = 0; j < numReachNodesEachSource; ++j) {
			v = nodesInIncreasingOrder[j];
			sourcesIn2OutIndex[i][j] = v;
			sourcesOut2InIndex[v] = j;
		}

		degrees.clear();
		for (j = 0; j < numReachNodesEachSource; ++j) 
			degrees.push_back(successorsAdjLists.size(nodesInIncreasingOrder[j]));
		sourcesSuccessorsAdjList[i]->allocateMemory(degrees);
		for (j = 0; j < numReachNodesEachSource; ++j)
		{
			v = nodesInIncreasingOrder[j];
			sz = successorsAdjLists.size(v);
			adj = &successorsAdjLists.access(v, 0);
			for (k = 0; k < sz; ++k)
				sourcesSuccessorsAdjList[i]->push_back(j, sourcesOut2InIndex[ adj[k] ] );
			sourcesNumShortestPathsAfter[i][j] = sourcesNumShortestPathsBefore[i][v];
			numEdgesInDAG += sourcesSuccessorsAdjList[i]->size(j);
		}
		for (j = 0; j < numReachNodesEachSource; ++j) // make them same
			sourcesNumShortestPathsBefore[i][j] = sourcesNumShortestPathsAfter[i][j];
	}	

	numEdgesInDAG /= L;
	clock_t ft1 = clock();
	printf("n: %d, numEdgesInDAG:%.0lf, ratio: %.2lf, time for finding shortest paths: %.3lf,\n", n, numEdgesInDAG, numEdgesInDAG  / n, (ft1 - st) / 1000.0);
	foundTimes.push_back((ft1 - st) / 1000.0);

	int innerLargestID;
	int *In2OutIndex;
	// find top-k nodes with the largest GBC from the L sources
	for (i = 0; i < n; ++i) isFoundforGBC[i] = false;

	double maxMemoryDis = 0;
	double sumGBC;
	double* gbcEachItself = new double[n];

	for (k = 0; k < K; ++k)
	{
		for (i = 0; i < n; ++i) betweenCentralities[i] = 0;  // outer indices
		for (i = 0; i < L; ++i)
		{
			src = sourceNodes[i];
			if (true == isFoundforGBC[src])
				continue; // this source node is already in the set of found nodes

			for (j = 0; j < n; ++j) currrentBetweenCentralities[j] = 0; // inner indices

			shortestPathsAdjLists = sourcesSuccessorsAdjList[i];
			oneSourceNumShortestPathsBefore = sourcesNumShortestPathsBefore[i];
			oneSourceNumShortestPathsAfter = sourcesNumShortestPathsAfter[i];
			for (j = 0; j < n; ++j) {
				if (0 == oneSourceNumShortestPathsAfter[j]) gbcEachItself[j] = 0;
				else gbcEachItself[j] = double(oneSourceNumShortestPathsAfter[j]) / oneSourceNumShortestPathsBefore[j];
			}	

			for (j = n - 1; j >= 0; --j)
			{
				if (0 == oneSourceNumShortestPathsAfter[j]) continue;

				sumGBC = 0;
				sz = shortestPathsAdjLists->size(j);
				adj = &shortestPathsAdjLists->access(j, 0);
				for (l = 0; l < sz; ++l)
				{
					w = adj[l];
					if (0 == oneSourceNumShortestPathsAfter[w]) continue;
					sumGBC += oneSourceNumShortestPathsAfter[j] * currrentBetweenCentralities[w] / oneSourceNumShortestPathsAfter[w];
				}
				//if (0 == k && 0 == i) maxMemoryDis += abs(adj[0] - j);
				currrentBetweenCentralities[j] = sumGBC + gbcEachItself[j];
				//currrentBetweenCentralities[j] +=  double(oneSourceNumShortestPathsAfter[j])/ oneSourceNumShortestPathsBefore[j];
			}
			currrentBetweenCentralities[0]--;
			

			In2OutIndex = sourcesIn2OutIndex[i];
			for (j = 0; j < n; ++j) betweenCentralities[ In2OutIndex[j] ] += currrentBetweenCentralities[j];
		}

		largestIncreasedBC = 0;
		outLargestNodeID = -1;
		for (j = 0; j < n; ++j)
			if (betweenCentralities[j] > largestIncreasedBC)
			{
				largestIncreasedBC = betweenCentralities[j];
				outLargestNodeID = j;
			}
		assert(outLargestNodeID >= 0);

		for (j = 0; j < n; ++j) decreaseDelta[j] = 0;
		for (j = 0; j < n; ++j) visited[j] = false;
		for (i = 0; i < L; ++i)
		{
			if (true == isFoundforGBC[sourceNodes[i]])
				continue; // this source node is already in the set of found nodes

			oneSourceNumShortestPathsAfter = sourcesNumShortestPathsAfter[i];
			shortestPathsAdjLists = sourcesSuccessorsAdjList[i];

			In2OutIndex = sourcesIn2OutIndex[i];
			for (j = 0; j < n; ++j)
				if (In2OutIndex[j] == outLargestNodeID) break;
			assert(j != n);
			innerLargestID = j;

			decreaseDelta[innerLargestID] = oneSourceNumShortestPathsAfter[innerLargestID];
			if (0 == decreaseDelta[innerLargestID] ) continue;

			visited[innerLargestID] = true;
			numRecordNodes = 1;
			recordNodes[0] = innerLargestID;

			queue_clear(); // clear the queue
			queue_push(innerLargestID); // the kth found node
			while (queue_size() > 0) // bfs, decrease the number of shorestest paths
			{
				u = queue_pop();
				sz = shortestPathsAdjLists->size(u);
				adj = &shortestPathsAdjLists->access(u, 0);
				for (j = 0; j < sz; ++j)
				{
					v = adj[j];
					if (0 == oneSourceNumShortestPathsAfter[v] ) continue;
					if (false == visited[v]) {
						queue_push(v);
						visited[v] = true;
						recordNodes[numRecordNodes] = v;
						numRecordNodes++;
					}
					decreaseDelta[v] += decreaseDelta[u];
				}
			}

			for (j = 0; j < numRecordNodes; ++j) {
				v = recordNodes[j];
				visited[v] = false;
				sourcesNumShortestPathsAfter[i][v] -= decreaseDelta[v];
				decreaseDelta[v] = 0;
			}
		}

		//if ((k + 1) % 10 == 0)
		printf("find the %d th node: %d, increased GBC: %.0lf, used total time: %.2lf, delta Time: %.3lf \n",	k + 1, outLargestNodeID, largestIncreasedBC, (clock() - st) / 1000.0, (clock() - ft1) / 1000.0);
		ft1 = clock();
		if ((k + 1) % 10 == 0) foundTimes.push_back((ft1 - st) / 1000.0);

		isFoundforGBC[outLargestNodeID] = true;
		foundNodes.push_back(outLargestNodeID);
		totalGBC += largestIncreasedBC;

		//	printf("%d th node: %d, increased BC: %.3lf\n", k + 1, 	largestNodeID, largestIncreasedBC);
	}// end for (int k = 0; k < K; ++k)

	maxMemoryDis /= L;
	maxMemoryDis /= n;
	printf("avergae dis: %.0lf\n", maxMemoryDis);
	totalGBC = (totalGBC / L);
	totalGBC /= double(n - 1.0); // normalized GBC

	for (i = 0; i < L; ++i) {
		delete sourcesSuccessorsAdjList[i];
		delete[] sourcesNumShortestPathsBefore[i];
		delete[]sourcesNumShortestPathsAfter[i];
		delete[]sourcesIn2OutIndex[i];
		
		
		sourcesSuccessorsAdjList[i] = NULL;
		sourcesNumShortestPathsBefore[i] = NULL;
		sourcesNumShortestPathsAfter[i] = NULL;
		sourcesIn2OutIndex[i] = NULL;
	}
	delete[]sourcesOut2InIndex;
	delete[]gbcEachItself;
	sourcesOut2InIndex = NULL;

	return totalGBC;
}

double Graph::rankByGroupBCWithSampleSourcesFastPlus(vector<int>& sourceNodes, vector<int>& foundNodes, int K, vector<double>& foundTimes)
{ // time complexity: L m + KL 1.5n
	assert(K > 0);
	foundNodes.clear();
	foundNodes.reserve(K);
	foundTimes.clear();
	int L = sourceNodes.size();

	int i, j, k, l;
	bool debug = false;
	nonUniformArray<int>  successorsAdjLists; // record successors in the shortest paths
	vector<int> degrees; // degree of each node 
	degrees.reserve(n);
	for (i = 0; i < n; ++i) degrees.push_back(adjLists.size(i));
	successorsAdjLists.allocateMemory(degrees);

	assert(L <= maxSampledSources && n <= MAX_ROWS_OneCC);
	for (i = 0; i < L; ++i) {
		sourcesSuccessorsAdjList[i] = new nonUniformArraySmaller<int>;
		sourcesNumShortestPathsBefore[i] = new int[n];
		sourcesNumShortestPathsAfter[i] = new int[n];
	}

	int decreaseDelta[MAX_ROWS];
	for (i = 0; i < L; ++i)
		for (j = 0; j < n; ++j)
			sourcesNumShortestPathsBefore[i][j] = 0;
	int recordNodes[MAX_ROWS + 1000];
	int numRecordNodes;

	int src;
	int u, v, w;
	// the number of nodes that their shortest distances to src have been found
	int sz;
	int* adj;

	double totalGBC = 0;
	double largestIncreasedBC;
	int largestNodeID = 0;

	//printf("num nodes: %d, numedges: %d\n", n, numOfEdges() *2 );
	clock_t st = clock();

	double numEdgesInDAG = 0;
	int newDistance;
	double maxMemoryDis = 0;
	for (i = 0; i < L; ++i) // find the shortest path DAG from each source
	{
		if (i % 10 == 0) printf("finding the %dth shortest-path DAG...\n", i + 1);
		src = sourceNodes[i];

		sourcesNumShortestPathsBefore[i][src] = 1;
		for (j = 0; j < n; ++j) successorsAdjLists.clear(j);
		// find the shortest paths from each node 
		for (j = 0; j < n; ++j) visited[j] = false;
		visited[src] = true;
		distances[src] = 0;
		queue_clear(); // clear the queue
		queue_push(src);
		while (queue_size() > 0) // bfs
		{
			u = queue_pop();

			sz = adjLists.size(u);
			adj = &adjLists.access(u, 0);
			newDistance = distances[u] + 1;
			for (j = 0; j < sz; ++j)
			{
				v = adj[j];
				if (false == visited[v])
				{
					visited[v] = true;
					distances[v] = newDistance;
					queue_push(v);
				}
				if (distances[v] == newDistance) {
					successorsAdjLists.push_back(u, v);
					sourcesNumShortestPathsBefore[i][v] += sourcesNumShortestPathsBefore[i][u];
				}
			}
		}
		degrees.clear();
		for (j = 0; j < n; ++j) degrees.push_back(successorsAdjLists.size(j));
		sourcesSuccessorsAdjList[i]->allocateMemory(degrees);
		for (j = 0; j < n; ++j)
		{
			sz = successorsAdjLists.size(j);
			adj = &successorsAdjLists.access(j, 0);
			for (k = 0; k < sz; ++k)
				sourcesSuccessorsAdjList[i]->push_back(j, adj[k]);
			numEdgesInDAG += sz;
			if (sz > 0) maxMemoryDis += abs(sourcesSuccessorsAdjList[i]->access(j,0)-j);
		}
	}
	for (i = 0; i < L; ++i)
		for (j = 0; j < n; ++j)
			sourcesNumShortestPathsAfter[i][j] = sourcesNumShortestPathsBefore[i][j];
	
	clock_t ft1 = clock();
	printf("n: %d, numEdgesInDAG:%.0lf, ratio: %.2lf, time for finding shortest paths: %.3lf,\n", n,  numEdgesInDAG / L, (numEdgesInDAG / L)/n, (ft1 -st)/1000.0);
	foundTimes.push_back( (ft1-st)/1000.0 );
	maxMemoryDis /= L;
	maxMemoryDis /= n;
	printf("maxMemoryDis: %.0lf\n", maxMemoryDis);

	//vector<int> nodesInIncreasingOrder;
	//nodesInIncreasingOrder.reserve(n);
	// find top-k nodes with the largest GBC from the L sources
	for (i = 0; i < n; ++i) isFoundforGBC[i] = false;
	for (k = 0; k < K; ++k)
	{
		//if (k % 5  == 0) {
			
		//}
		for (i = 0; i < n; ++i) betweenCentralities[i] = 0;
		for (i = 0; i < L; ++i)
		{
			src = sourceNodes[i];
			if (true == isFoundforGBC[src])
				continue; // this source node is already in the set of found nodes

			for (j = 0; j < n; ++j) currrentBetweenCentralities[j] = 0;
			for (j = 0; j < n; ++j) visited[j] = false;
			visited[src] = true;
			for (j = 0; j < foundNodes.size(); ++j) visited[foundNodes[j]] = true;

			shortestPathsAdjLists = sourcesSuccessorsAdjList[i];
			oneSourceNumShortestPathsBefore = sourcesNumShortestPathsBefore[i];
			oneSourceNumShortestPathsAfter = sourcesNumShortestPathsAfter[i];
			/*queue_clear(); // clear the queue
			queue_push(src);
			nodesInIncreasingOrder.clear();
			nodesInIncreasingOrder.push_back(src);
			while (queue_size() > 0) // bfs
			{
				u = queue_pop();
				sz = shortestPathsAdjLists->size(u);
				adj = &shortestPathsAdjLists->access(u, 0);
				for (j = 0; j < sz; ++j)
				{
					v = adj[j];
					if (true == isFoundforGBC[v]) continue;
					if (false == visited[v])
					{
						visited[v] = true;
						queue_push(v);
						nodesInIncreasingOrder.push_back(v);
					}
				}
			}
			for (j = nodesInIncreasingOrder.size() - 1; j >= 0; --j)
			{
				v = nodesInIncreasingOrder[j];
				sz = shortestPathsAdjLists->size(v);
				adj = &shortestPathsAdjLists->access(v, 0);
				for (l = 0; l < sz; ++l)
				{
					w = adj[l];
					if (true == isFoundforGBC[w]) continue;
					currrentBetweenCentralities[v] += oneSourceNumShortestPathsAfter[v]* currrentBetweenCentralities[w] / oneSourceNumShortestPathsAfter[w];
				}
				currrentBetweenCentralities[v] += (double)oneSourceNumShortestPathsAfter[v] / oneSourceNumShortestPathsBefore[v];
			}/**/
			recursiveCalGBC(src);
			currrentBetweenCentralities[src]--;
			for (j = 0; j < n; ++j)
				betweenCentralities[j] += currrentBetweenCentralities[j];
		}
		largestIncreasedBC = 0;
		largestNodeID = -1;
		//largestNodeID = k+50;
		for (j = 0; j < n; ++j)
			if (betweenCentralities[j] > largestIncreasedBC)
			{
				largestIncreasedBC = betweenCentralities[j];
				largestNodeID = j;
			}
		assert(largestNodeID >= 0);
		
		for (j = 0; j < n; ++j) decreaseDelta[j] = 0;
		for (j = 0; j < n; ++j) visited[j] = false;
		for (i = 0; i < L; ++i)
		{
			if (true == isFoundforGBC[ sourceNodes[i] ])
				continue; // this source node is already in the set of found nodes
			
			visited[largestNodeID] = true;
			numRecordNodes = 1;
			recordNodes[0] = largestNodeID;

			for (j = 0; j < foundNodes.size(); ++j) 
				visited[foundNodes[j]] = true;
			
			decreaseDelta[ largestNodeID ] = sourcesNumShortestPathsAfter[i][largestNodeID];
			queue_clear(); // clear the queue
			queue_push( largestNodeID ); // the kth found node
			while (queue_size() > 0) // bfs, decrease the number of shorestest paths
			{
				u = queue_pop();
				sz = sourcesSuccessorsAdjList[i]->size(u);
				adj = &sourcesSuccessorsAdjList[i]->access(u, 0);
				for (j = 0; j < sz; ++j)
				{
					v = adj[j];
					if (true == isFoundforGBC[v]) continue;
					if (false == visited[v]) {
						queue_push(v);
						visited[v] = true;
						recordNodes[numRecordNodes] = v;
						numRecordNodes++;
					}
					decreaseDelta[v] += decreaseDelta[u];
				}
			}
			
			for (j = 0; j < numRecordNodes; ++j) {
				v = recordNodes[j];
				visited[v] = false;
				sourcesNumShortestPathsAfter[i][v] -= decreaseDelta[v];
				decreaseDelta[v] = 0;
				//if (0 == sourcesNumShortestPathsAfter[i][v] ) sourcesSuccessorsAdjList[i].clear(v);
			}
			for (j = 0; j < foundNodes.size(); ++j)
				visited[foundNodes[j]] = false;
		}

		//if ((k + 1) % 10 == 0) 
			printf("find the %d th node: %d, increased GBC: %.2lf, used total time: %.2lf, delta Time: %.2lf \n", 	k + 1, largestNodeID, largestIncreasedBC, (clock() - st) / 1000.0, (clock() - ft1) / 1000.0);
		ft1 = clock();
		if( (k+1)%10== 0 ) foundTimes.push_back((ft1 - st) / 1000.0);

		isFoundforGBC[largestNodeID] = true;
		foundNodes.push_back(largestNodeID);
		totalGBC += largestIncreasedBC;
		//	printf("%d th node: %d, increased BC: %.3lf\n", k + 1, 	largestNodeID, largestIncreasedBC);
	}// end for (int k = 0; k < K; ++k)
	//clock_t ft2 = clock();
	//printf("\ntime used, phase 1: %d ms, phase 2: %d ms\n\n", (int)(ft1 - st), (int)(ft2 - ft1));

	totalGBC = (totalGBC / L);
	totalGBC /= double(n - 1.0); // normalized GBC
	
	for (i = 0; i < L; ++i) {
		delete sourcesSuccessorsAdjList[i];
		delete[] sourcesNumShortestPathsBefore[i];
		delete []sourcesNumShortestPathsAfter[i];
		sourcesSuccessorsAdjList[i] = NULL;
		sourcesNumShortestPathsBefore[i] = NULL;
		sourcesNumShortestPathsAfter[i] = NULL;
	}
	
	return totalGBC;
}


/*double Graph::rankByGroupBCWithSampleSourcesFast(vector<int>& sourceNodes, vector<int>& foundNodes, int K)
{ // time complexity: L m + KL 2n
	assert(K > 0);
	foundNodes.clear();
	foundNodes.reserve(K);
	int L = sourceNodes.size();

	int i, j, k, l;
	bool debug = false;
	nonUniformArray<int>  successorsAdjLists; // record successors in the shortest paths
	vector<int> degrees; // degree of each node 
	degrees.reserve(n);
	for (i = 0; i < n; ++i) degrees.push_back(adjLists.size(i));
	successorsAdjLists.allocateMemory(degrees);

	//const int maxSampledSources = 300;
	//nonUniformArray<int> sourcesSuccessorsAdjList[maxSampledSources];
	assert(L <= maxSampledSources && n <= MAX_ROWS_OneCC);
	//int  sourcesNumShortestPathsBefore[maxSampledSources][MAX_ROWS];
	for (i = 0; i < L; ++i)
		for (j = 0; j < n; ++j)
			sourcesNumShortestPathsBefore[i][j] = 0;

	int src;
	int u, v, w;
	// the number of nodes that their shortest distances to src have been found
	int sz;
	int* adj;

	double totalGBC = 0;
	double largestIncreasedBC;
	int largestNodeID = 0;

	vector<int> nodesInIncreasingOrder;
	nodesInIncreasingOrder.reserve(n);

	//printf("num nodes: %d, numedges: %d\n", n, numOfEdges() *2 );
	double numEdgesInDAG = 0;
	int newDistance;
	for (i = 0; i < L; ++i) // find the shortest path DAG from each source
	{
		src = sourceNodes[i];

		sourcesNumShortestPathsBefore[i][src] = 1;
		for (j = 0; j < n; ++j) successorsAdjLists.clear(j);
		// find the shortest paths from each node 
		for (j = 0; j < n; ++j) visited[j] = false;
		visited[src] = true;
		distances[src] = 0;
		queue_clear(); // clear the queue
		queue_push(src);
		while (queue_size() > 0) // bfs
		{
			u = queue_pop();

			sz = adjLists.size(u);
			adj = &adjLists.access(u, 0);
			newDistance = distances[u] + 1;
			for (j = 0; j < sz; ++j)
			{
				v = adj[j];
				if (false == visited[v])
				{
					visited[v] = true;
					distances[v] = newDistance;
					queue_push(v);
				}
				if (distances[v] == newDistance) {
					successorsAdjLists.push_back(u, v);
					sourcesNumShortestPathsBefore[i][v] += sourcesNumShortestPathsBefore[i][u];
				}
			}
		}
		degrees.clear();
		for (j = 0; j < n; ++j) degrees.push_back(successorsAdjLists.size(j));
		sourcesSuccessorsAdjList[i].allocateMemory(degrees);
		for (j = 0; j < n; ++j)
		{
			sz = successorsAdjLists.size(j);
			adj = &successorsAdjLists.access(j, 0);
			for (k = 0; k < sz; ++k)
				sourcesSuccessorsAdjList[i].push_back(j, adj[k]);
			numEdgesInDAG += sz;
		}
	}
	numEdgesInDAG /= L;
	//printf("n: %d, numEdgesInDAG:%.0lf, times: %.2lf\n", n, numEdgesInDAG, numEdgesInDAG / n);

	// find top-k nodes with the largest GBC from the L sources
	for (i = 0; i < n; ++i) isFoundforGBC[i] = false;
	for (k = 0; k < K; ++k)
	{
		for (i = 0; i < n; ++i) betweenCentralities[i] = 0;
		for (i = 0; i < L; ++i)
		{
			src = sourceNodes[i];
			if (true == isFoundforGBC[src])
				continue; // this source node is already in the set of found nodes
			
			for (j = 0; j < n; ++j) numShortestPathsAfter[j] = 0;
			numShortestPathsAfter[src] = 1;
			queue_clear(); // clear the queue
			queue_push(src);
			nodesInIncreasingOrder.clear();
			nodesInIncreasingOrder.push_back(src);
			
			for (j = 0; j < n; ++j) visited[j] = false;
			visited[src] = true;
			for (j = 0; j < foundNodes.size(); ++j) visited[foundNodes[j]] = true;
			while (queue_size() > 0) // bfs
			{
				u = queue_pop();
				sz = sourcesSuccessorsAdjList[i].size(u);
				adj = &sourcesSuccessorsAdjList[i].access(u, 0);

				for (j = 0; j < sz; ++j)
				{
					v = adj[j];
					if (false == visited[v])
					{
						visited[v] = true;
						queue_push(v);
						numShortestPathsAfter[v] += numShortestPathsAfter[u];
						nodesInIncreasingOrder.push_back(v);
					}
					else if (false == isFoundforGBC[v])
						numShortestPathsAfter[v] += numShortestPathsAfter[u];
				}
			}
			
			for (j = 0; j < n; ++j) currrentBetweenCentralities[j] = 0;
			for (j = nodesInIncreasingOrder.size() - 1; j >= 0; --j)
			{
				v = nodesInIncreasingOrder[j];
				sz = sourcesSuccessorsAdjList[i].size(v);
				adj = &sourcesSuccessorsAdjList[i].access(v, 0);
				for (l = 0; l < sz; ++l)
				{
					w = adj[l];
					if (true == isFoundforGBC[w]) continue;
					currrentBetweenCentralities[v] += numShortestPathsAfter[v]
						* currrentBetweenCentralities[w] / numShortestPathsAfter[w];
				}
				currrentBetweenCentralities[v] += (double)numShortestPathsAfter[v] / sourcesNumShortestPathsBefore[i][v];
			}
			currrentBetweenCentralities[src]--;

			
			for (j = 0; j < n; ++j)
				betweenCentralities[j] += currrentBetweenCentralities[j];
		}// end for(i = 0; i < L; ++i)
		largestIncreasedBC = 0;
		for (j = 0; j < n; ++j)
		{
			//if(betweenCentralities[j] > 0) assert(isFoundforGBC[j] == false);
			if (betweenCentralities[j] > largestIncreasedBC)
			{
				largestIncreasedBC = betweenCentralities[j];
				largestNodeID = j;
			}
		}

		//assert(largestIncreasedBC > 0);
		isFoundforGBC[largestNodeID] = true;
		foundNodes.push_back(largestNodeID);
		totalGBC += largestIncreasedBC;
		//	printf("%d th node: %d, increased BC: %.3lf\n", k + 1, 	largestNodeID, largestIncreasedBC);
	}// end for (int k = 0; k < K; ++k)

	totalGBC = (totalGBC / L) * n;
	return totalGBC / (n * (n - 1));
}/**/


double Graph::rankByGroupBCWithSampleSources(vector<int>& sourceNodes, vector<int>& foundNodes, int K)
{ // time complexity: O(K L m)
	foundNodes.clear();
	assert(K > 0);
	int L = sourceNodes.size();

	int i, j, k, l;
	bool debug = false;
	nonUniformArray<int>  successorsAdjLists; // record successors in the shortest paths
	vector<int> degrees; // degree of each node 
	degrees.reserve(n);
	for (i = 0; i < n; ++i) degrees.push_back(adjLists.size(i));
	successorsAdjLists.allocateMemory(degrees);
	int src;
	int u, v, w;
	// the number of nodes that their shortest distances to src have been found
	int sz;
	int* adj;

	double totalGBC = 0;
	double largestIncreasedBC;
	int largestNodeID=0;

	vector<int> nodesInIncreasingOrder;
	nodesInIncreasingOrder.reserve(n);

	//printf("num nodes: %d, numedges: %d\n", n, numOfEdges() *2 );
	int newDistance;
	for (i = 0; i < n; ++i) isFoundforGBC[i] = false;
	clock_t st = clock();
	for (k = 0; k < K; ++k)
	{
		//printf("finding the %dth node....\n", k + 1);
		for (i = 0; i < n; ++i) betweenCentralities[i] = 0;
		for (i = 0; i < L; ++i)
		{
			src = sourceNodes[i];
			if (true == isFoundforGBC[src] ) 
				continue; // this source node is already in the set of found nodes

			for (j = 0; j < n; ++j) successorsAdjLists.clear(j);
			for (j = 0; j < n; ++j) numShortestPathsBefore[j] = 0;
			numShortestPathsBefore[src] = 1;
			// find the shortest paths from each node 
			for (j = 0; j < n; ++j) visited[j] = false;
			visited[src] = true;
			distances[src] = 0;
			queue_clear(); // clear the queue
			queue_push(src);
			while (queue_size() > 0) // bfs
			{
				u = queue_pop();

				sz = adjLists.size(u);
				adj = &adjLists.access(u, 0);
				newDistance = distances[u] + 1;
				for (j = 0; j < sz; ++j)
				{
					v = adj[j];
					if (false == visited[v])
					{
						visited[v] = true;
						distances[v] = newDistance;
						queue_push(v);
					}
					if (distances[v] == newDistance) {
						numShortestPathsBefore[v] += numShortestPathsBefore[u];
						successorsAdjLists.push_back(u, v);
					}
				}
			}
			/*if (src == 0 && debug)
			{
				printf("shortest path DAG from source %d:----------------\n", src);
				successorsAdjLists.printArray();
				printf("numShortestPathsBefore:\n");
				for (j = 0; j < n; ++j) printf("%d,", numShortestPathsBefore[j] );
				printf("\n");
			}	/**/

			for (j = 0; j < n; ++j) numShortestPathsAfter[j] = 0;
			numShortestPathsAfter[src] = 1;
			queue_clear(); // clear the queue
			queue_push(src);
			nodesInIncreasingOrder.clear();
			nodesInIncreasingOrder.push_back(src);
			for (j = 0; j < n; ++j) visited[j] = false;
			visited[src] = true;
			for (j = 0; j < foundNodes.size(); ++j) visited[ foundNodes[j] ] = true;
			while (queue_size() > 0) // bfs
			{
				u = queue_pop();
				sz = successorsAdjLists.size(u);
				adj = &successorsAdjLists.access(u, 0);
				
				for (j = 0; j < sz; ++j)
				{
					v = adj[j];
					if (false == visited[v])
					{
						visited[v] = true;
						queue_push(v);
						numShortestPathsAfter[v] += numShortestPathsAfter[u];
						nodesInIncreasingOrder.push_back(v);
					}else if (false == isFoundforGBC[v])
						numShortestPathsAfter[v] += numShortestPathsAfter[u];
				}
			}
			/*if (src == 0 && debug)
			{
				for (j = 0; j < n; ++j)
					printf("to node %d, numShortestPathsRes: %d\n", j, numShortestPathsAfter[j]);
			}/**/
			for (j = 0; j < n; ++j) currrentBetweenCentralities[j] = 0;
			for (j = nodesInIncreasingOrder.size() - 1; j >= 0; --j)
			{ 
				v = nodesInIncreasingOrder[j];
				sz = successorsAdjLists.size(v);
				adj = &successorsAdjLists.access(v, 0);
				for (l = 0; l < sz; ++l)
				{
					w = adj[l];
					if (true == isFoundforGBC[w]) continue;
					currrentBetweenCentralities[v] += numShortestPathsAfter[v]* currrentBetweenCentralities[w] / numShortestPathsAfter[w];
				}
				currrentBetweenCentralities[v] +=(double)numShortestPathsAfter[v] / numShortestPathsBefore[v];
			}
			currrentBetweenCentralities[src]--; 

			/*if (src == 0 && debug)
			{
				for (j = 0; j < n; ++j)
					printf("from source %d, node %d, BC: %.3lf\n", src, j, 
						currrentBetweenCentralities[j]);
			}/**/
			for (j = 0; j < n; ++j)
				betweenCentralities[j] += currrentBetweenCentralities[j];
			//if ((i + 1) % 50 == 0) printf("finding the %d th node, %d th source, used time: %.3lf\n", k+1, i+1, (clock()-st)/1000.0);
		}// end for(i = 0; i < L; ++i)
		largestIncreasedBC = 0;
		for (j = 0; j < n; ++j)
		{
			//if(betweenCentralities[j] > 0) assert(isFoundforGBC[j] == false);
			if (betweenCentralities[j] > largestIncreasedBC)
			{
				largestIncreasedBC = betweenCentralities[j];
				largestNodeID = j;
			}
		}
			
		//assert(largestIncreasedBC > 0);
		isFoundforGBC[largestNodeID] = true;
		foundNodes.push_back(largestNodeID);
		totalGBC += largestIncreasedBC;
	//	printf("%d th node: %d, increased BC: %.3lf\n", k + 1, 	largestNodeID, largestIncreasedBC);
	}// end for (int k = 0; k < K; ++k)

	totalGBC /= L;
	totalGBC /= double(n - 1.0);
	return totalGBC;
}



void Graph::rankByBetweennessCentrality(vector<int>& rankedNodes)
{ //Pre: the graph is connected and the weight of each node is one
	// Return: nodes in the network are ranked by their BC scores
	// usage: Brandes' Alglorithm
	clock_t st = clock();
	int i, j, k;
	bool debug = false;
	nonUniformArray<int>  successorsAdjLists; // record successors in the shortest paths
	vector<int> degrees; // degree of each node 
	degrees.reserve(n);
	for (i = 0; i < n; ++i) degrees.push_back( adjLists.size(i) );
	successorsAdjLists.allocateMemory(degrees);

	int src;
	int u, v, w;
	// the number of nodes that their shortest distances to src have been found
	//int numDistancesFound = 0;
	int sz;
	int* adj;

	double* betweenCentralities = new double[n];
	for (i = 0; i < n; ++i) betweenCentralities[i] = 0;
	double* currrentBetweenCentralities = new double[n];

	vector<int> nodesInIncreasingOrder;
	nodesInIncreasingOrder.reserve(n);

	int newDistance;
	int initialEdges = 0;
	int totalEdges = 0;
	
	for (i = 0; i < n; ++i)
	{
		//if (i % 100 == 0 && i > 0) printf("%d th node\n", i);
		src = i;
		//if (src == 0) debug = true;
		//else		  debug = false;

		nodesInIncreasingOrder.clear();
		nodesInIncreasingOrder.push_back(src);

		for (j = 0; j < n; ++j) successorsAdjLists.clear(j);

		for (j = 0; j < n; ++j) numShortestPathsBefore[j] = 0;
		numShortestPathsBefore[src] = 1;

		 // find the shortest paths from each node 
		for (j = 0; j < n; ++j) visited[j] = false;
		visited[src] = true;
		distances[src] = 0;
		queue_clear(); // clear the queue
		queue_push(src);

		while (queue_size() > 0)
		{
			u = queue_pop();

			sz = adjLists.size(u);
			adj = &adjLists.access(u, 0);
			newDistance = distances[u] + 1;
			for (j = 0; j < sz; ++j)
			{
				v = adj[j];
				initialEdges++;
				if (false == visited[v])
				{
					visited[v] = true;
					distances[v] = newDistance;
					queue_push(v);
					nodesInIncreasingOrder.push_back(v);
					
				}
				if (distances[v] == newDistance) {
					numShortestPathsBefore[v] += numShortestPathsBefore[u];
					successorsAdjLists.push_back(u, v);
				}
			}
		}
		/*if (debug) {
			printf("Number of shortest paths:\n");
			for (j = 0; j < n; ++j)
				printf("%d, numShortestPaths: %d\n", j, numShortestPaths[j]);
			printf("\nSuccessors of each node in the shortest paths:\n");
			successors.printGraphTopology();
		}/**/

		for (j = 0; j < n; ++j) currrentBetweenCentralities[j] = 0;
		for (j = nodesInIncreasingOrder.size() - 1; j > 0; --j)
		{// the source node is excluded, which is at the 1st location of the array
			v = nodesInIncreasingOrder[j];
			sz = successorsAdjLists.size(v);
			adj = &successorsAdjLists.access(v, 0);
			for (k = 0; k < sz; ++k)
			{
				w = adj[k];
				currrentBetweenCentralities[v] += numShortestPathsBefore[v]
					* (1.0 + currrentBetweenCentralities[w]) / numShortestPathsBefore[w];
			}
		}/**/

		/*if (debug) {
			printf("\nBetweenness centrality of each node for source node: %d\n", src);
			for (j = 0; j < n; ++j)
				printf("node: %d, BC: %.3lf\n", j,
					currrentBetweenCentralities[j]);
		}/**/
		for (j = 0; j < n; ++j)
			betweenCentralities[j] += currrentBetweenCentralities[j];
	}

	printf("n: %d, initialEdges: %d, totalEdges: %d\n", n, initialEdges / n, totalEdges / n);

	vector<Node> allNodes;
	allNodes.reserve(n);
	for (i = 0; i < n; ++i)
		allNodes.push_back(Node(i, betweenCentralities[i]));
	sort(allNodes.begin(), allNodes.end(), comp); //sort the count in decreasing order
	Node node;
	rankedNodes.clear();
	rankedNodes.reserve(n);
	for (i = 0; i < allNodes.size(); ++i)
	{
		node = allNodes[i];
		rankedNodes.push_back(node.id);
		//if (i < 10) printf("%d, BC: %.0lf\n", i+1, node.key);
	}

	delete[]betweenCentralities;
	delete[]currrentBetweenCentralities;
}


void Graph::shortestDistances(int src, bool markSPT)
{
	int i, j;
	for (i = 0; i < n; ++i) visited[i] = false;

	visited[src] = true;
	distances[src] = 0;
	queue_clear(); // clear the queue
	queue_push(src);

	int u, v;
	// the number of nodes that their shortest distances to src have been found
	//int numDistancesFound = 0;
	int sz;
	int *adj; 
	while (queue_size() > 0 )
	{
		u = queue_pop();

		//++numDistancesFound;
		//if (debug) printf("node %d, dis: %d\n", u, (int)distances[u]);

		sz = adjLists.size(u);
		adj = & adjLists.access(u, 0);
		for (j = 0; j < sz; ++j)
		{
			//v = adjLists.access(u, j);
			if ( visited[ adj[j] ] ) continue;
			v = adj[ j ];
			visited[v] = true;
			distances[v] = distances[u] + 1;
			queue_push(v);
		}
	}
}

void Graph::allocateVariablesForDFS()
{
	assert(n > 0);
	if (NULL == ccBelongto)
	{
		numChildren = new int[n];
		numNodesSubtree = new int[n];
		parent = new int[n];
		isAP = new bool[n];
		discoveredTime = new int[n];
		lowestTime = new int[n];
		ccBelongto = new int[n];
		lowerBounds = new double[n];
	}

	for (int i = 0; i < n; ++i) // initialize
	{
		visited[i] = false;
		numChildren[i] = 0;
		numNodesSubtree[i] = 0;
		parent[i] = -1; // -1 indicates that the node is the root
		isAP[i] = false;
		lowerBounds[i] = 0;
	}
	time = 0;
}


void Graph::initilizeVariables()
{
	numChildren = NULL;
	numNodesSubtree = NULL;
	parent = NULL;
	isAP = NULL;
	discoveredTime = NULL;
	lowestTime = NULL;
	ccBelongto = NULL;
	lowerBounds = NULL;
	time = 0;
}

void Graph::clearVariables()
{
	if (NULL != ccBelongto)
	{
		delete []numChildren;
		delete []numNodesSubtree;
		delete []parent;
		delete []isAP;
		delete []discoveredTime;
		delete []lowestTime;
		delete []ccBelongto;
		delete []lowerBounds;
	}
	initilizeVariables();
}

void Graph::findAPsAndLowerBounds()
{
	nonUniformArray<int> allCCs;
	findAllCCsDFS(allCCs);
	
	int numNodesValid = 0; 
	int i;
	for (i = 0; i < allCCs.size(); ++i)
		numNodesValid += allCCs.size(i);

	allocateVariablesForDFS();

	for (i = 0; i < n; ++i)
	{
		if (false == isValid[i]) continue;
		if (true == visited[i]) continue;
		
		modifiedDFS(i, allCCs, numNodesValid);
		//if ( sizeOfCC > 10 ) printf("CC %d, size: %d\n", numofCCs, sizeOfCC);
	}
	
	double base = 0;
	for (i = 0; i < allCCs.size(); ++i)
		base += allCCs.size(i) * (numNodesValid - allCCs.size(i));
	base -= 2 * numNodesValid;

	//printf("the set of APs:\n");
	int numOfAPs = 0;
	for (i = 0; i < n; ++i)
	{
		if (false == isValid[i])continue;
		//if (false == isAP[i]) continue;
		++numOfAPs;
		lowerBounds[i] += base; 

		lowerBounds[i] *= zeta;
		//printf("%d ", i);
	}
	//printf("\n");
	//printf("num of APs: %d, percentage: %.2lf \n",  numOfAPs, 100.0*numOfAPs / n);

	/*printf("\nDFS tree: \n");
	for (i = 0; i < n; ++i)
		printf("parent of %d: %d\n", i, parent[i]);

	printf("\n");
	for (i = 0; i < n; ++i)
		printf("node %d, numofchildren: %d, numofdesc: %d\n", i, numChildren[i], numNodesSubtree[i]);

	printf("\n");
	for (i = 0; i < n; ++i)
		printf("node %d, discoverTime: %d, lowestTime: %d\n", i, discoveredTime[i], lowestTime[i]);
    */
}

int Graph::findAllCCsDFS(nonUniformArray<int> &allCCs)
{
	int i; 
	int numOfCCs = 0;
	sizeofEachCC.clear();
	sizeofEachCC.reserve( 1000 );

	allocateVariablesForDFS();

	int sizeOfCC;
	int maxCCsize = 0;
	for (i = 0; i < n; ++i)
	{
		if (false == isValid[i]) continue;
		if (true == visited[i]) continue;

		sizeOfCC = 0;
		DFS(i, numOfCCs, sizeOfCC); // dfs search
		//assert(sizeOfCC > 0);
		sizeofEachCC.push_back( sizeOfCC );
		++numOfCCs;
		//if ( sizeOfCC >= 30 )   printf("CC %d, size: %d\n", numOfCCs, sizeOfCC);
		if ( sizeOfCC > maxCCsize) maxCCsize = sizeOfCC;
	}
	if (debug) 
		printf("numofCCs: %d, maxCCsize: %d (%.2lf)\n", numOfCCs, maxCCsize, 100.0 * maxCCsize / n);
	
	if (debug)
		for (i = 0; i < n; ++i) printf("node: %d, in CC %d\n", i, ccBelongto[i]);

	allCCs.allocateMemory( sizeofEachCC );
	for (i = 0; i < n; ++i)
	{
		if (false == isValid[i]) continue; 
		allCCs.push_back(ccBelongto[i], i);
	}

	if (debug)
		for (i = 0; i < allCCs.size(); ++i)
		{
			printf("CC %d: ", i);
			for (int j = 0; j < allCCs.size(i); ++j)
				printf("%d ", allCCs.access(i, j ));
			printf("\n");
		}
	return numOfCCs;
}

int Graph::findAllCCsBFS(nonUniformArray<int> &allCCs)
// find all connected components
{
	int i, j;
	int numofCCs = 0;
	sizeofEachCC.clear();
	sizeofEachCC.reserve(1000);

	allocateVariablesForDFS();

	queue_clear();

	int sizeOfCC;
	int maxCCsize = 0;
	int u, v;
	int sz;
	int *adj;
	for (i = 0; i < n; ++i)
	{
		if (true == visited[i]) continue;
		if (false == isValid[i]) continue;

		// search the cc that contains vertex i
		sizeOfCC = 1;
		ccBelongto[ i ] = numofCCs;
		queue_push(i);
		while (queue_size() > 0 )
		{
			u = queue_pop();

			sz = adjLists.size(u);
			adj = & adjLists.access(u, 0);
			for (j = 0; j < sz; ++j)
			{
				v = adj[ j ];
				if (true == visited[v]) continue;
				if (false == isValid[v]) continue;

				queue_push(v);

				visited[ v ] = true;
				ccBelongto[ v ] = numofCCs;
				++sizeOfCC;
			}
		}
		assert(sizeOfCC > 0);

		sizeofEachCC.push_back(sizeOfCC);
		++numofCCs;
		//if ( sizeOfCC > 10 )   printf("CC %d, size: %d\n", numofCCs, sizeOfCC);
		if (sizeOfCC > maxCCsize) maxCCsize = sizeOfCC;
	}
	if (debug)
		printf("numofCCs: %d, maxCCsize: %d (%.2lf)\n", numofCCs, maxCCsize, 100.0 * maxCCsize / n);

	if (debug)
		for (i = 0; i < n; ++i) printf("node: %d, in CC %d\n", i, ccBelongto[i]);

	allCCs.allocateMemory( sizeofEachCC );
	for (i = 0; i < n; ++i)
	{
		if (false == isValid[i]) continue;
		allCCs.push_back(ccBelongto[i], i);
	}

	if (debug)
		for (i = 0; i < allCCs.size(); ++i)
		{
		printf("CC %d: ", i);
		for (int j = 0; j < allCCs.size(i); ++j)
			printf("%d ", allCCs.access(i, j));
		printf("\n");
		}
	return numofCCs;
}


inline void Graph::DFS(int u, int CCID, int & sizeOfCC)
{
	visited[u] = true;
	ccBelongto[u] = CCID;
	++sizeOfCC;

	int v;
	int sz = adjLists.size(u);
	int *adj = & adjLists.access(u, 0);
	for (int j = 0; j < sz ; ++j)
	{
		v = adj[ j ];
		//assert(v >= 0 && v < n);

		if (isValid[v] && !visited[v] )
		{
			parent[v] = u;
			DFS(v, CCID, sizeOfCC);
			numNodesSubtree[u] += numNodesSubtree[v];
		}
	}
	++numNodesSubtree[u]; // including itself
}

void Graph::modifiedDFS(int u, nonUniformArray<int> & allCCs, int numNodesValid)
{
	visited[ u ] = true;
	++time;
	discoveredTime[ u ] = lowestTime[ u ] = time;
	//printf("time %d, node %d\n", time, u);
	int ccSize = allCCs.size(ccBelongto[u]);
	int cc0 = ccSize - 1;

	int v;
	int j;
	int sz = adjLists.size(u);
	int *adj = &adjLists.access(u, 0);
	int nv; 
	for (j = 0; j < sz; ++j)
	{
		v = adj[j];
		if (false == isValid[v]) continue;

		if (false == visited[v] )
		{
			parent[ v ] = u;
			++numChildren[ u ];

			modifiedDFS(v, allCCs, numNodesValid);
			numNodesSubtree[u] += numNodesSubtree[v];
			
			if (lowestTime[v] < lowestTime[u])
				lowestTime[u] = lowestTime[v];
			// u is an articulation point in following cases

			// (1) u is root of DFS tree and has two or more chilren.
			if (parent[u] == -1 && numChildren[u] >= 2 ||
				// (2) If u is not root and low value of one of its child
				// is more than discovery value of u.
				parent[u] != -1 && lowestTime[v] >= discoveredTime[u])
			{
				isAP[u] = true;

				nv = numNodesSubtree[v];
				assert(nv < ccSize);
				lowerBounds[u] += nv * ( ccSize - 1 - nv);
				cc0 -= nv;
				assert(cc0 >= 0);
			}
		}
		else if (v != parent[u]) // a back edge is found
		{
			if (discoveredTime[v] < lowestTime[u])
				lowestTime[u] = discoveredTime[v];
		}
	}
	++numNodesSubtree[u]; // including itself

	if (true == isAP[u])
	{
		lowerBounds[u] += cc0 * (ccSize - 1 - cc0); // *zeta;
		lowerBounds[u] += 2 * ccSize;
	}else  // that is the upper bound for non-APs
		lowerBounds[u] = 2 * ccSize + 1; // *zeta;
	
}

void Graph::printGraph()
{
	int i, j;
	for (i = 0; i < n; ++i)
	{
		printf("node %d, weight:%f,  adjList: \n", i, nodesWeights[i]);
		for (j = 0; j < adjLists.size(i); ++j)
			printf("(%d, %d, %.4lf)\t", i, adjLists.access(i,j), adjWeights.access(i,j) );
		printf("\n\n");
	}
}

void Graph::printGraphTopology()
{
	int i, j;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < adjLists.size(i); ++j)
			printf("(%d, %d)\t", i, adjLists.access(i, j));
		printf("\n\n");
	}
}

void Graph::cleanGraph(const char *graphFile)
{
	int i, j;
	ifstream in(graphFile);
	if (!in) {
		printf("Graph file `%s' was not found.\n", graphFile);
		return;
	}
	clearGraph();

	in >> n >> m;
	assert(n > 0 && m >= 0);
	assert(n <= MAX_NODES);

	vector<int> degrees; // the degree of each node 
	degrees.resize(n, 0);

	vector<Edge> allEdges;
	allEdges.reserve(m);

	int u, v;
	int numEdges = 0;
	for (i = 0; i < m; ++i)
	{
		if (i % 1000000 == 0) printf("%d (1000,000) th  edge \n", (int)(i / 1000000));

		in >> u >> v;

		if (u == v) {
			printf("A self loop occurs.\n");  
			continue;
		}
		if (u > v) swap(u, v); // make sure that u < v

		if (u < 0 || u >= n){ printf("u exceeds range: %d\n", u); continue; }
		if (v < 0 || v >= n){ printf("v exceeds range: %d\n", v); continue; }

		++degrees[u];
		allEdges.push_back(Edge(u, v));
	}
	in.close();

	adjLists.allocateMemory( degrees );
	for (i = 0; i < allEdges.size(); ++i)
	{
		u = allEdges[i].u;
		v = allEdges[i].v;
		assert(u < v);
		if (true == adjLists.find(u, v)) continue;
		//printf("Parallel edge (%d, %d) occurs\n", u, v); 
		adjLists.push_back(u, v); // insert vertex v to the adjlist of vertex u
		++numEdges;
	}
	printf("number of edges (before): %d, (after): %d\n", m, numEdges);
	m = numEdges;

	char dstFile[100];
	sprintf_s(dstFile, "datasets/%s", graphFile);
	ofstream out(dstFile);
	if (!out) abort();

	out << n << '\t' << m << '\n';
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < adjLists.size(i); ++j)
			out << i << '\t' << adjLists.access(i,j) << '\n';
	}
	out.close();

	printf("clean finished\n");
}


void Graph::generateEdgeNodeWeightsCaseStudy(double minNodeWeight, double maxNodeWeight)
{
	assert(0 <= minNodeWeight && minNodeWeight <= maxNodeWeight);
	vector<int> degrees; // degree of each node
	degrees.reserve(numOfNodes());
	int i, j;
	for (i = 0; i < adjLists.size(); ++i)
		degrees.push_back(adjLists.size(i));
	adjWeights.allocateMemory(degrees);
	double w;
	int u, v;
	for (i = 0; i < adjLists.size(); ++i)
	{
		u = i;
		for (j = 0; j < adjLists.size(i); ++j)
		{
			v = adjLists.access(i, j);
			if (u == 9 || v == 9)	w = 0.1;
			else w = 0.2;
			
			adjWeights.push_back(i, w);
			//if (i <= 1) printf("(%d,%d, %f)\n", i, adjLists.access(i, j), w);
		}
	}


	nodesWeights.clear();
	nodesWeights.reserve(numOfNodes());
	for (j = 0; j < numOfNodes(); ++j)
	{
		w = Uniform() *(maxNodeWeight - minNodeWeight) + minNodeWeight;
		nodesWeights.push_back(w);
		//if (j < 10) printf("(%d, %f)\n", j, w);
	}
	//printGraph();
}

void Graph::generateEdgeNodeWeights(double minEdgeWeight, double maxEdgeWeight, double minNodeWeight, double maxNodeWeight)
{
	assert(0 <= minEdgeWeight&& minEdgeWeight <= maxEdgeWeight && maxEdgeWeight <= 1.0);
	assert(0 <= minNodeWeight && minNodeWeight <= maxNodeWeight);
	vector<int> degrees; // degree of each node
	degrees.reserve(numOfNodes());
	int i, j;
	for (i = 0; i < adjLists.size(); ++i)
		degrees.push_back(adjLists.size(i));
	adjWeights.allocateMemory(degrees);
	double w;
	int degree;
	for (i = 0; i < adjLists.size(); ++i)
	{
		degree = adjLists.size(i);
		for (j = 0; j < adjLists.size(i); ++j)
		{
			w = Uniform() * (maxEdgeWeight - minEdgeWeight) + minEdgeWeight;
			adjWeights.push_back(i, w);
			//if (i <= 1) printf("(%d,%d, %f)\n", i, adjLists.access(i, j), w);
		}
	}
	
	nodesWeights.clear();
	nodesWeights.reserve(numOfNodes());
	for (j = 0; j < numOfNodes(); ++j)
	{
		w = Uniform() *(maxNodeWeight - minNodeWeight) + minNodeWeight;
		nodesWeights.push_back(w);
		//if (j < 10) printf("(%d, %f)\n", j, w);
	}
}

void Graph::readGraph(const char *graphFile, bool isDirected)
{
	FILE *f; 
	errno_t err = fopen_s(&f, graphFile, "r");
	if (err != 0)
	{
		printf("file: %s was not found!\n", graphFile);
		abort();
	}

	clearGraph();

	fscanf_s(f, "%d%d", &n, &m);
	
	assert(n > 0 && m >= 0);
	assert(n <= MAX_NODES);

	isValid = new bool[n];
	int i;
	for (i = 0; i < n; ++i) isValid[i] = true;

	vector<Edge> allEdges;
	allEdges.reserve(m);

	vector<int> degrees; // degree of each node
	degrees.resize(n, 0);

	int u, v;
	int numEdges = 0;
	for (i = 0; i < m; ++i)
	{
		if (i % 1000000 == 0 && i > 0 ) printf("%d (1000,000) th  edge \n", (int)(i / 1000000));

		fscanf_s(f, "%d%d", &u, &v);
		assert(u >= 0 && u < n);
		assert(v >= 0 && v < n);

		
		++numEdges;

		++degrees[ u ];
		if (false == isDirected) ++degrees[v];
		allEdges.push_back( Edge(u, v) );
	}
	fclose(f);
	assert(m == numEdges);

	adjLists.allocateMemory(degrees);
	assert( n == adjLists.size() );

	for (i = 0; i < allEdges.size(); ++i)
	{
		u = allEdges[i].u; 
		v = allEdges[i].v;
		adjLists.push_back(u, v);
		if (false == isDirected)
			adjLists.push_back(v, u);
	}
	vector<int> adjListSorting;
	adjListSorting.reserve(n);
	int j;
	for (i = 0; i < n; ++i)
	{
		adjListSorting.clear();
		for (j = 0; j < adjLists.size(i); ++j)
			adjListSorting.push_back(adjLists.access(i, j));
		sort(adjListSorting.begin(), adjListSorting.end());// sort in ascending order
		adjLists.clear(i);
		for (j = 0; j < adjListSorting.size(); ++j)
			adjLists.push_back(i, adjListSorting[j]);
	}
	printf("graph has been created.\n");
}

void Graph::clearGraph()
{
	if (NULL != isValid)  delete[]isValid;
	isValid = NULL;
	adjLists.clear();
	n = m = 0;
	clearVariables();
}

Graph::Graph()
{
	isValid = NULL;
	n = m = 0;
	debug = false;

	initilizeVariables();
}

Graph::~Graph()
{
	clearGraph();
	clearVariables();
}

Graph::Graph(const Graph &g)
{
	assert(g.n > 0);
	n = g.n; 
	m = g.m;
	adjLists = g.adjLists;
	isValid = new bool[n];
	for (int i = 0; i < n; ++i) isValid[i] = g.isValid[i];
	debug = g.debug;
	
	initilizeVariables();
}

Graph & Graph::operator=(const Graph &g)
{
	if (this == &g) return *this;

	clearGraph();
	clearVariables();

	assert(g.n > 0);
	n = g.n;
	m = g.m;
	adjLists = g.adjLists;

	isValid = new bool[n];
	for (int i = 0; i < n; ++i) isValid[i] = g.isValid[i];
	debug = g.debug;

	initilizeVariables();
	
	return *this;
}



