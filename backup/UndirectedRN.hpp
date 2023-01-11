#pragma once


// *** INCLUDES *** //
#include <iostream>
#include <vector>
#include <set>
#include <unordered_map>

#include "GeographicLib/Geodesic.hpp"

#include "lemon/list_graph.h"
#include "lemon/path.h"



class UndirectedRN
{
public:

	// *** PUBLIC TYPEDEFS *** //

	typedef typename lemon::ListGraph Graph;	// Type of the graph used internally.
	typedef typename std::pair<double,double> Spatial2DPoint; // Container for a geopoint (lat,long).
	typedef typename std::pair<double, 			// Type used to represent POIs (i.e., cost...
							   uint32_t> POI;	// + ID graph node.



protected:

	// *** PROTECTED CLASSES *** //

	// Generic comparator used for priority queues.
	template<typename K>
	struct CmpPQueue {inline bool operator()(const K& lhs, const K& rhs) const {return lhs.first > rhs.first;}};



	// *** PROTECTED FIELDS *** //

	Graph g;

	std::vector<Graph::Node> vecNodes;
	Graph::NodeMap<Spatial2DPoint> mapVertices;

	std::vector<uint32_t> startPosTypes; std::vector<uint32_t> endPosTypes;
	std::vector<uint32_t> startPosTypesBackup; std::vector<uint32_t> endPosTypesBackup;

	std::vector<std::unordered_map<uint32_t, double>> vecCostPOIs;
	std::vector<std::unordered_map<uint32_t, double>> vecCostPOIsBackup;

	std::vector<std::vector<POI>> sortedVecCostPOIs;
	std::vector<std::vector<POI>> sortedVecCostPOIsBackup;

	std::vector<Graph::Edge> vecEdges;
	Graph::EdgeMap<double> mapEdges;

	// Fields used by the baseline variant employing precomputed network distances
	// between POIs belonging to adjacent COIs.
	std::unordered_map<uint32_t, std::vector<std::pair<double, uint32_t>>> tableTDs;
	std::unordered_map<uint32_t, std::vector<std::vector<std::pair<double, uint32_t>>>> tableTypeTDs;



	// *** PROTECTED METHODS *** //

	void parseFileVertices(const char* namefile);
	void parseFileEdges(const char* namefile);

	void parseFilePOIs(const std::vector<const char*>& vecNameFilesPOIs);
	void modifyEdges(const std::vector<std::set<std::pair<double,Graph::Node>>>& vecModEdges);

	inline double computeDistanceVertices(const Graph::Node& p1, const Graph::Node& p2);



	// *** PROTECTED CTORS *** //

	UndirectedRN() : mapVertices(g), mapEdges(g) {};



public:

	// *** PUBLIC CTORS/DTOR *** //

	UndirectedRN(const char* nameFileVertices,
	 	 	 	 const char* nameFileEdges,
				 std::vector<const char*>& vecNameFilesPOIs);
	~UndirectedRN() {};



	// *** PUBLIC METHODS *** //

	// Methods dealing with the computation of shortest paths.
	inline static double computeDistancePoints(const Spatial2DPoint& p1, const Spatial2DPoint& p2);
	inline double computeDistancePointsIDs(const uint32_t& p1, const uint32_t& p2);
	double SPDijkstra(const uint32_t& s, const uint32_t& d);
	double SPAStar(const uint32_t& s, const uint32_t& d);

	// Methods implementing INE and PNE.
	std::pair<int32_t, double> INELookup(const uint32_t& idSource, const uint32_t& typePOI, const uint32_t& k = 0);
	std::pair<std::vector<uint32_t>, double> PNE(const uint32_t& idSource);

	inline uint32_t getNumVertices() const {return this->vecNodes.size();};
	inline uint32_t getNumEdges() const {return this->vecEdges.size();};

	// Methods dealing with POI categories.
	inline double getCostPOI(const uint32_t& type, const uint32_t& idPOI) {return(this->vecCostPOIs[type][idPOI]);};
	inline double getCostSequence(const std::vector<uint32_t>& seq);
	inline uint32_t getNumTypesPOIs() {return this->sortedVecCostPOIs.size();};
	inline const std::vector<std::vector<POI>>& getVecSortedPOIs() {return this->sortedVecCostPOIs;};

	void setupQuerySequence(const uint32_t& lenSequence);

	// Methods used to precompute Euclidean and network distances bewteen POIs.
	void precomputeEDPOIs(const std::string& outputFile);
	void precomputeTDPOIs(const std::string& outputFile);

	// Methods used to load precomputed network distances between POIs.
	void loadPrecomputedTDPOIs(const std::string& nameFile);

	// Methods concerning the baseline using precomputed network distances between POIs.
	std::pair<int32_t, double> INEPrecomp(const uint32_t& idSource, const uint32_t& typePOI, const uint32_t& k = 0);
	void PNEBaseline(const uint32_t& idSource, const double& minCost);
	void resetHashTableTDs() {this->tableTypeTDs.clear();};
	void PNEBaselinePrecomp(const uint32_t& idSource, const double& minCost);
	void PNEBaselinePrecompConvSkyline(const uint32_t& idSource, const double& minCost);
	void PNEBaselinePrecompWithoutPruning(const uint32_t& idSource, const double& minCost);

	void initPrecomputedNNsPOIs();
	void initkNNsSource(const uint32_t& idSource);

	// Test methods.
	void testPNE();
};



// *** PROTECTED METHODS DEFINITIONS *** //

double UndirectedRN::computeDistanceVertices(const Graph::Node& p1, const Graph::Node& p2)
{
	return(UndirectedRN::computeDistancePoints(this->mapVertices[p1], this->mapVertices[p2]));
}



// *** PUBLIC METHODS DEFINITIONS *** //

/**
 * @brief This method computes the distance between two geographical points expressed as a pair (lat,lon).
 * @note The distance returned is expressed in meters.
 */
double UndirectedRN::computeDistancePoints(const Spatial2DPoint& p1, const Spatial2DPoint& p2)
{
	static const GeographicLib::Geodesic& geod(GeographicLib::Geodesic::WGS84());
	double distance_metres = 0;
	geod.Inverse(p1.first, p1.second,
				 p2.first, p2.second, distance_metres);

	return(distance_metres);
}

/**
 * @brief This method computes the Euclidean distance between two points.
 */
double UndirectedRN::computeDistancePointsIDs(const uint32_t& p1, const uint32_t& p2)
{
	return(UndirectedRN::computeDistancePoints(this->mapVertices[this->g.nodeFromId(p1)],
											   this->mapVertices[this->g.nodeFromId(p2)]));
}

double UndirectedRN::getCostSequence(const std::vector<uint32_t>& seq)
{
	double cost = 0;
	for(uint32_t type = 0; type < seq.size(); type++)
		cost += this->vecCostPOIs[type][seq[type]];
	return(cost);
}
