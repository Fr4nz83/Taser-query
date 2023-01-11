// *** INCLUDES *** //

#include "URN.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>



// *** PROTECTED METHODS DEFINITIONS *** //

/**
 * @brief This method copies a graph (i.e., the set of data structures that represent it) into another one.
 */
void URN::copyGraph(Graph& gOld, std::vector<Graph::Node>& vecNodesOld, Graph::NodeMap<Spatial2DPoint>& mapVerticesOld, std::vector<Graph::Edge>& vecEdgesOld, Graph::EdgeMap<double>& mapEdgesOld,
					Graph& gNew, std::vector<Graph::Node>& vecNodesNew, Graph::NodeMap<Spatial2DPoint>& mapVerticesNew, std::vector<Graph::Edge>& vecEdgesNew, Graph::EdgeMap<double>& mapEdgesNew)
{
	// Instantiante copy object.
	lemon::GraphCopy<Graph, Graph> cg(gOld, gNew);

	// Copy node map
	cg.nodeMap(mapVerticesOld, mapVerticesNew);

	// Copy edge map
	cg.edgeMap(mapEdgesOld, mapEdgesNew);

	// Copy nodes
	vecNodesNew = std::vector<Graph::Node>(vecNodesOld.size());
	for(uint32_t i = 0; i < vecNodesOld.size(); i++)
		cg.node(vecNodesOld[i], vecNodesNew[i]);

	// Copy edges
	vecEdgesNew = std::vector<Graph::Edge>(vecEdgesOld.size());
	for(uint32_t i = 0; i < vecEdgesOld.size(); i++)
		cg.edge(vecEdgesOld[i], vecEdgesNew[i]);

	cg.run();
}

/**
 * @brief This method parses the POIs from the files specified in "vecNameFiles" and stores their information in an appropriate vector.
 * 		  The graph is NOT modified at this stage!
 */
void URN::parseFilePOI(const std::vector<const char*>& vecNameFiles)
{
	std::cout << "INIT => Number of POI files to load: " << vecNameFiles.size() << "\n";

	std::ifstream infile;
	uint32_t countPOI = this->vecNodes.size();
	for(auto namefile : vecNameFiles)
	{
		infile.open(namefile, std::ios::binary);
		if(!infile.is_open())
		{
			std::cout << "Error opening file: \"" << namefile << "\". Terminating...\n";
			exit(1);
		}


		std::cout << "INIT => Reading POI file: " << namefile << "\n";
		this->vecCOIsBackup.push_back({});
		uint32_t idEdge; double distStart; double cost;
		while(infile >> cost >> cost >> cost >> idEdge >> distStart >> cost)
		{
			if(idEdge < this->vecEdges.size())
			{
				this->vecCOIsBackup.back().push_back({idEdge, distStart, cost, countPOI});
				countPOI++;
			}
			// Verify that there aren't POIs that refer to non-existing edges.
			else
			{
				std::cout << "Error parsing POI! Terminating!\n";
				exit(0);
			}
		}


		infile.close();
		std::cout << "INIT => Number of POIs belonging to type " << this->vecCOIsBackup.size() - 1 << ": " << this->vecCOIsBackup.back().size() << "\n";
	}

	std::cout << "INIT => Number of COIs: " << this->vecCOIsBackup.size() << "\n";
}

/**
 * @brief Randomly extract up to "numPOIsCOI" from each available COI. Subsequently, it updates the graph and initializes all the appropriate data structures.
 *
 * @param numPOIsCOI Maximum number of POIs allowed in a COI.
 */
void URN::applyPOIs(uint32_t numPOIsCOI)
{
	if(numPOIsCOI == 0) numPOIsCOI = std::numeric_limits<uint32_t>::max();
	std::cout << "Maximum number of POIs per COI: " << numPOIsCOI << "\n";


	// 1 - Extract random POIs from COIs.
	this->vecCOIs.clear();
	for(auto tmp : this->vecCOIsBackup)
	{
		std::random_shuffle(tmp.begin(), tmp.end());
		if(tmp.size() > numPOIsCOI) tmp.resize(numPOIsCOI);
		this->vecCOIs.push_back(tmp);
	}


	// 2 - Create the mapping between old IDs and new IDs of the POIs in the TD hash table.
	std::unordered_map<uint32_t, uint32_t> tmpMap;
	uint32_t newIndexPOI = this->vecNodes.size();
	for(const auto& vecPOI : this->vecCOIs)
		for(const auto& p : vecPOI)
			tmpMap[p.idOrgVertex] = newIndexPOI++;


	// 3 - Reinitialize TD hash table with the new mapping.
	this->tableTDs.clear();
	for(const auto& tmp : tmpMap)
	{
		this->tableTDs[tmp.second] = std::vector<std::pair<double, uint32_t>>();
		for(auto old : this->tableTDsBackup[tmp.first])
		{
			// If the POI in the original hash table is present among the selected ones, proceed to add it to the new
			// hash table.
			if(tmpMap.find(old.second) != tmpMap.end())
			{
				old.second = tmpMap[old.second];
				this->tableTDs[tmp.second].push_back(old);
			}
		}
	}


	// 4 - Update the graph with the set of POIs obtained above.
	this->updateGraph();

}

/**
 * @brief This method updates the graph with COIs having their POIs extracted randomly from the original COIs.
 */
void URN::updateGraph()
{
	// Reset main data structures graph.
	this->startPosTypes.clear();
	this->endPosTypes.clear();
	this->vecCostPOIs.clear();
	this->sortedVecCostPOIs.clear();


	std::cout << "INIT => Number of COIs: " << this->vecCOIs.size() << "\n";
	std::vector<std::set<std::pair<double,Graph::Node>>> vecModEdges; vecModEdges.resize(this->vecEdges.size());
	for(const auto& vecPOI : this->vecCOIs)
	{
		uint32_t countPOI = 0;
		std::cout << "INIT => Starting position POI type " << this->startPosTypes.size() << ": " << this->vecNodes.size() << " in the vertices' vector\n";

		std::set<POI> tmpRes;
		this->vecCostPOIs.push_back(std::unordered_map<uint32_t, double>());
		this->startPosTypes.push_back(this->vecNodes.size()); // Save the starting position associated with this type of POIs.
		for(const auto& p : vecPOI)
		{
			if(p.idEdge < vecModEdges.size())
			{
				const uint32_t indexCOI = this->startPosTypes.size() - 1; // Get the index of the current COI.
				const uint32_t indexPOI = this->vecNodes.size();		  // Get the index of the vertex associated with this POI.
				this->vecCostPOIs[indexCOI][indexPOI] = p.cost;			  // Store the POI's cost in the unordered map.
				tmpRes.insert(POI(p.cost, indexPOI));						  // Store the POI's cost in an ordered set.

				// Add POI to the graph and save the object representing it in vecNodes.
				this->vecNodes.push_back(g.addNode());

				// std::cout << "Edge to modify: " << idEdge << "\n";
				// Take note of the POI to be associated with this edge;
				vecModEdges[p.idEdge].insert(std::make_pair(p.distStart, this->vecNodes.back()));
				countPOI++;
			}
			// NOTE: there are some edgeIDs in POI files that do not actually exist!
			//       Take care of these occurrences.
			else
			{
				std::cout << "Error parsing POI!\n";
			}
		}
		// Take note of the last position of the POI belonging to this COI.
		this->endPosTypes.push_back(this->vecNodes.size() - 1);
		// Sort the POIs of the current COI according to their costs and add them to the vector of sorted POIs.
		this->sortedVecCostPOIs.push_back(std::vector<POI>(tmpRes.begin(), tmpRes.end()));

		std::cout << "INIT => Number of POIs belonging to type " << this->startPosTypes.size() - 1 << ": " << countPOI << "\n";
	}
	// *** Update the edges in which we have at least one POI *** //
	this->modifyEdges(vecModEdges);

	// Save in backup variables the original content of POIs.
	this->vecCostPOIsBackup = this->vecCostPOIs;
	this->sortedVecCostPOIsBackup = this->sortedVecCostPOIs;
	this->startPosTypesBackup = this->startPosTypes;
	this->endPosTypesBackup = this->endPosTypes;
}



// *** PUBLIC CTORS *** //

/**
 * @brief This static method parse the content of a road network from a set of textual files
 * 		  to create an in-memory graph-based representation.
 */
URN::URN(const char* nameFileVertices,
		 const char* nameFileEdges,
		 std::vector<const char*>& vecNameFilesPOIs) :

UndirectedRN(),
mapVerticesBackup(gBackup),
mapEdgesBackup(gBackup)
{
	// Parse the file containing the vertices.
	this->parseFileVertices(nameFileVertices);

	// Parse the file containing the edges.
	this->parseFileEdges(nameFileEdges);

    // Copy the untouched road network into a set of backup variables.
	this->copyGraph(this->g, this->vecNodes, this->mapVertices, this->vecEdges, this->mapEdges,
					this->gBackup, this->vecNodesBackup, this->mapVerticesBackup, this->vecEdgesBackup, this->mapEdgesBackup);

	// Parse the file containing the POIs.
	this->parseFilePOI(vecNameFilesPOIs);


	// DEBUG.
	std::cout << "INIT => Number of vertices original road network: " << this->vecNodes.size() << "\n";
	std::cout << "INIT => Number of edges original road network: " << this->vecEdges.size() << "\n";
}



// *** PUBLIC METHODS *** //

void URN::loadPrecomputedTDPOIs(const std::string& nameFile)
{
	std::ifstream infile;
	infile.open(nameFile, std::ios::binary);
	if(!infile.is_open())
	{
		std::cout << "Error opening TDs file: \"" << nameFile << "\". Terminating...\n";
		exit(1);
	}

	// Init the hash tables.
	uint32_t source, target;
	double dist;
	while(infile >> source >> target >> dist)
	{
		if(this->tableTDsBackup.find(source) == this->tableTDsBackup.end())
			this->tableTDsBackup[source] = std::vector<std::pair<double, uint32_t>>();
		this->tableTDsBackup[source].push_back(std::make_pair(dist, target));

		if(this->tableTDsBackup.find(target) == this->tableTDsBackup.end())
			this->tableTDsBackup[target] = std::vector<std::pair<double, uint32_t>>();
		this->tableTDsBackup[target].push_back(std::make_pair(dist, source));
	}
	infile.close();

	// Sort the POIs associated with each POI in increasing order of travel distance.
	for(auto& poi : this->tableTDsBackup)
		std::sort(poi.second.begin(), poi.second.end());
}

/**
 * @brief This method resets the graph to the original road network and then applies COIs for which the POIs are extracted randomly
 *        from an initial set of COIs.
 */
void URN::regenerateGraph(const uint32_t& numPOIsCOI)
{
	// Restore the original graph.
	this->copyGraph(this->gBackup, this->vecNodesBackup, this->mapVerticesBackup, this->vecEdgesBackup, this->mapEdgesBackup,
					this->g, this->vecNodes, this->mapVertices, this->vecEdges, this->mapEdges);

	// Extract randomly some POIs and add them to the graph.
	this->applyPOIs(numPOIsCOI);
}

