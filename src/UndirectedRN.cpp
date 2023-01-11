// *** INCLUDES *** //

#include "UndirectedRN.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <random>
#include <queue>
#include <chrono>
#include "omp.h"

#include <boost/heap/fibonacci_heap.hpp>

#include "GeographicLib/GeodesicLine.hpp"

#include "SimpleSkyline.hpp"
#include "SimpleConventionalSkyline.hpp"
#include "lemon/dijkstra.h"
#include "lemon/astar.h"
// #include "lemon/abistar.h"

#include "EstimateCalculator.hpp"



// *** PROTECTED METHODS DEFINITIONS *** //

/**
 * @brief This method takes in input a file containing the vertices of a road network and
 * 		  store them in appropriate data structures in memory.
 */
void UndirectedRN::parseFileVertices(const char* namefile)
{
	// Open the input game trace file;
	std::ifstream infile; infile.open(namefile, std::ios::binary);

	if(!infile.is_open())
	{
		std::cout << "Error opening file: \"" << namefile << "\". Terminating...\n";
		exit(1);
	}


	// Initialize from file the list of vertices of the graph.
	std::cout << "INIT => Reading vertices file: \"" << namefile << "\"\n";
	Spatial2DPoint tmp;
	while(infile >> tmp.first >> tmp.first >> tmp.second)
	{
		// Add node to the graph and save the object representing the vertex.
		this->vecNodes.push_back(g.addNode());

		// Save lat, lon of the vertex in the map.
		this->mapVertices[this->vecNodes.back()] = tmp;
	}

	infile.close();
	std::cout << "INIT => Number of vertices: " << this->vecNodes.size() << "\n";
}

/**
 * @brief This method takes in input a file containing the edges of a road network and
 * 		  store them in appropriate data structures in memory.
 */
void UndirectedRN::parseFileEdges(const char* namefile)
{
	// Open the input game trace file;
	std::ifstream infile; infile.open(namefile, std::ios::binary);

	if(!infile.is_open())
	{
		std::cout << "Error opening file: \"" << namefile << "\". Terminating...\n";
		exit(1);
	}


	// Initialize from file the list of vertices of the graph.
	std::cout << "INIT => Reading edges file: \"" << namefile << "\"\n";
	uint32_t v1, v2;
	while(infile >> v1 >> v1 >> v2)
	{
		// Add edge to the graph;
		this->vecEdges.push_back(g.addEdge(this->vecNodes[v1], this->vecNodes[v2]));

		// Associate weight to the edge.
		this->mapEdges[this->vecEdges.back()] = this->computeDistancePoints(this->mapVertices[this->vecNodes[v1]], this->mapVertices[this->vecNodes[v2]]);
	}


	infile.close();
	std::cout << "INIT => Number of edges: " << this->vecEdges.size() << "\n";
}

/**
 * @brief This method parses the POIs from the files specified in "vecNameFiles" and update the
 * 		  set of vertices and the set of edges accordingly.
 */
void UndirectedRN::parseFilePOIs(const std::vector<const char*>& vecNameFiles)
{
	std::cout << "INIT => Number of POI files to load: " << vecNameFiles.size() << "\n";

	std::vector<std::set<std::pair<double,Graph::Node>>> vecModEdges; vecModEdges.resize(this->vecEdges.size());
	std::ifstream infile;
	for(auto namefile : vecNameFiles)
	{
		uint32_t countPOI = 0;

		infile.open(namefile, std::ios::binary);
		if(!infile.is_open())
		{
			std::cout << "Error opening file: \"" << namefile << "\". Terminating...\n";
			exit(1);
		}

		std::cout << "INIT => Reading POI file: " << namefile << "\n";
		std::cout << "INIT => Starting position POI type " << this->startPosTypes.size() << ": " << this->vecNodes.size() << " in the vertices' vector\n";

		std::set<POI> tmpRes;
		this->vecCostPOIs.push_back(std::unordered_map<uint32_t, double>());
		uint32_t idEdge; double distStart; double cost;
		this->startPosTypes.push_back(this->vecNodes.size()); // Save the starting position associated with this type of POIs.
		while(infile >> cost >> cost >> cost >> idEdge >> distStart >> cost)
		{
			if(idEdge < vecModEdges.size())
			{
				const uint32_t indexCOI = this->startPosTypes.size() - 1; // Get the index of the current COI.
				const uint32_t indexPOI = this->vecNodes.size();		  // Get the index of the vertex associated with this POI.
				this->vecCostPOIs[indexCOI][indexPOI] = cost;			  // Store the POI's cost in the unordered map.
				tmpRes.insert(POI(cost, indexPOI));						  // Store the POI's cost in an ordered set.

				// Add POI to the graph and save the object representing it in vecNodes.
				this->vecNodes.push_back(g.addNode());

				// std::cout << "Edge to modify: " << idEdge << "\n";
				// Take note of the POI to be associated with this edge;
				vecModEdges[idEdge].insert(std::make_pair(distStart, this->vecNodes.back()));
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


		infile.close();
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

/**
 * @brief This method makes the appropriate changes over the edges that contain at least one POI.
 */
void UndirectedRN::modifyEdges(const std::vector<std::set<std::pair<double,Graph::Node>>>& vecModEdges)
{
	static const GeographicLib::Geodesic& geod(GeographicLib::Geodesic::WGS84());
	for(uint32_t i = 0; i < vecModEdges.size(); i++)
	{
		// If there are modifications to do on this edge, proceed...
		if(vecModEdges[i].size())
		{
//			// DEBUG.
//			if(vecModEdges[i].size() >= 3)
//			{
//				std::cout << "Edge " << i << " (" << this->g.id(this->g.u(this->vecEdges[i]))
//						  << "," << this->g.id(this->g.v(this->vecEdges[i])) << ")\n";
//				std::cout << "Number of POIs on this edge: " << vecModEdges[i].size() << "\n";
//			}

			// *** ADD THE POIs TO THE EDGE *** //
			// (i) Retrieve the starting and ending nodes of the edge.
			Graph::Node startNode = this->g.u(this->vecEdges[i]);
			Graph::Node endNode = this->g.v(this->vecEdges[i]);
			// (ii) Retrieve their geo-coords...
			const Spatial2DPoint startEdge = this->mapVertices[startNode];
			const Spatial2DPoint endEdge = this->mapVertices[endNode];
			// (iii) Compute the line between them...
			GeographicLib::GeodesicLine line = geod.InverseLine(startEdge.first, startEdge.second,
																endEdge.first, endEdge.second);
			// ...and (iv) the distance separating them.
			const double lengthLine = line.Distance();


			// DEBUG.
//			if(vecModEdges[i].size() >= 3)
//				std::cout << "Length old edge: " << lengthLine << " metres \n";


			double lat = 0, lon = 0;
			auto it = vecModEdges[i].cbegin(); // Set up the iterator referring the POIs that have to
											   // be associated with the current edge.

			// *** Modify the original edge *** //
			this->g.changeV(this->vecEdges[i], it->second); // Update the ending vertex of the original edge.
			line.Position(it->first * lengthLine, lat, lon); // Compute the estimated coordinates of the POI.
			this->mapVertices[it->second] = Spatial2DPoint(lat,lon); // Update the coordinates of the vertex associated with the POI.
			this->mapEdges[this->vecEdges[i]] = it->first * lengthLine; // Update the weight of the edge.


			// DEBUG.
//			double cnt = 0;
//			if(vecModEdges[i].size() >= 3)
//			{
//				std::cout << "Length first edge: " << this->mapEdges[this->vecEdges[i]] << " metres \n";
//				cnt += this->mapEdges[this->vecEdges[i]];
//			}


			// *** Build a chain of edges that ends in the target vertex of the "old" edge. *** //
			Graph::Node lastNode = it->second;
			it++;
			while(it != vecModEdges[i].cend())
			{
				// Add new edge to the graph;
				this->vecEdges.push_back(g.addEdge(lastNode, it->second));

				// Set (i) the geocoords of the POI and (ii) the weight of the edge.
				line.Position(it->first * lengthLine, lat, lon); // Compute the estimated coordinates of the POI.
				this->mapVertices[it->second] = Spatial2DPoint(lat,lon); // Update the coordinates of the vertex associated with the POI.
				this->mapEdges[this->vecEdges.back()] = this->computeDistancePoints(this->mapVertices[lastNode], this->mapVertices[it->second]);


				// DEBUG.
//				if(vecModEdges[i].size() >= 3)
//				{
//					std::cout << "Length edge: " << this->mapEdges[this->vecEdges.back()] << " metres \n";
//					cnt += this->mapEdges[this->vecEdges.back()];
//				}


				// Update the node used to mark the end of the last edge.
				// Update the iterator.
				lastNode = it->second;
				it++;
			}

			// *** Add last edge of the chain (it must point to the end vertex of the "original" edge) *** //
			this->vecEdges.push_back(g.addEdge(lastNode, endNode));
			// Set weight of the edge.
			this->mapEdges[this->vecEdges.back()] = this->computeDistancePoints(this->mapVertices[lastNode], this->mapVertices[endNode]);


			// DEBUG.
//			if(vecModEdges[i].size() >= 3)
//			{
//				std::cout << "Length final edge: " << this->mapEdges[this->vecEdges.back()] << " metres \n";
//				cnt += this->mapEdges[this->vecEdges.back()];
//				std::cout << "Sum lengths: " << cnt << " metres \n";
//			}
		}
	}
}

/**
 * @brief This method implements the Progressive Network Expansion to compute the optimal sequenced route.
 *
 * @param idSource ID of the starting vertex.
 */
std::pair<std::vector<uint32_t>, double> UndirectedRN::PNE(const uint32_t& idSource)
{
	// Some typedefs and declarations.
	typedef struct
	{   std::vector<uint32_t> elSeq; std::vector<uint32_t> elK; std::vector<double> elW;
	} Sequence;
	typedef std::pair<double, Sequence> ElQueue;
	std::priority_queue<ElQueue, std::vector<ElQueue>, CmpPQueue<ElQueue>> queuePaths;


	// Explore the search space.
	const uint32_t numTypes = this->sortedVecCostPOIs.size();
	queuePaths.push(ElQueue(0, {std::vector<uint32_t>({idSource}), std::vector<uint32_t>({0}), std::vector<double>({0})}));
	ElQueue el;
	while(!queuePaths.empty())
	{
		// Pop the top element of the queue.
		el = queuePaths.top(); queuePaths.pop();

		// DEBUG.
//		std::cout << "Popped sequence ";
//		for(auto elem : el.second.elSeq)
//			std::cout << elem << " ";
//		std::cout << "(W " << el.first << ")\n";

		// 1 - Check if we have found the solution: if so, exit the cycle and return the sequence.
		if(el.second.elSeq.size() == numTypes + 1)
			break;


		// 2 - Otherwise, find out the NN with respect to the last COI of the partial sequence.
		const auto NN = this->INELookup(el.second.elSeq.back(),		  // Source node.
								  	    el.second.elSeq.size() - 1,   // Target category.
										el.second.elK.back());		  // Find k-th NN.


		// 3 - Generate the new partial sequence;
		Sequence tmpSeq = el.second;
		tmpSeq.elSeq.push_back(NN.first); tmpSeq.elW.push_back(tmpSeq.elW.back() + NN.second); tmpSeq.elK.push_back(0);
		queuePaths.push(ElQueue(el.first + NN.second, tmpSeq));

		// DEBUG.
//		std::cout << "Generated sequence: ";
//		for(auto elem : tmpSeq.elSeq)
//			std::cout << elem << " ";
//		std::cout << "(W " << el.first + NN.second << ")\n";


		// 4 - Look for the next NN for the predecessor of el. Then, add it to the queue.
		if(el.second.elSeq.size() > 1)
		{
			Sequence& seq = el.second;

			seq.elSeq.pop_back(); seq.elK.pop_back(); seq.elW.pop_back(); // Remove the last category from the sequence.
			seq.elK.back()++; 		   // Update the k of the next NN to find.
			el.first = seq.elW.back(); // Update the partial weight in el accordingly.

			// Find the k-th NN (if possible).
//			const auto kappa = seq.elK.back();
//			std::cout << "TEST: " << seq.elK.size() - 1 << "," << this->sortedVecCostPOIs[seq.elK.size() - 1].size() << "\n";
			if(seq.elK.back() < this->sortedVecCostPOIs[seq.elK.size() - 1].size())
			{
				const auto NNPred = this->INELookup(seq.elSeq.back(),		// Source node.
											  	    seq.elSeq.size() - 1,   // Target category.
													seq.elK.back());		// Find k-th NN.

				// Update el and push it in the queue.
				seq.elSeq.push_back(NNPred.first); seq.elW.push_back(seq.elW.back() + NNPred.second); seq.elK.push_back(0);
				el.first = seq.elW.back();
				queuePaths.push(el);

				// DEBUG.
//				std::cout << "Generated prec. sequence (k=" << kappa << "): ";
//				for(auto elem : el.second.elSeq)
//					std::cout << elem << " ";
//				std::cout << "(W " << el.first << ")\n";
			}
//			else
//				std::cout << "Reached END!!\n";
		}
	}

	// Return the OSR.
	return(std::pair<std::vector<uint32_t>, double>(el.second.elSeq, el.first));
}



// *** PUBLIC CTORS *** //

/**
 * @brief This static method parse the content of a road network from a set of textual files
 * 		  to create an in-memory graph-based representation.
 */
UndirectedRN::UndirectedRN(const char* nameFileVertices,
						   const char* nameFileEdges,
						   std::vector<const char*>& vecNameFilesPOIs) :

mapVertices(g),
mapEdges(g)
{

	// Parse the file containing the vertices.
	this->parseFileVertices(nameFileVertices);
	this->parseFileEdges(nameFileEdges);
	this->parseFilePOIs(vecNameFilesPOIs);


	// DEBUG.
	std::cout << "INIT => Number of vertices (after POIs): " << this->vecNodes.size() << "\n";
	std::cout << "INIT => Number of edges (after POIs): " << this->vecEdges.size() << "\n";
}



// *** PUBLIC METHODS *** //

/**
 * @brief This method generates a random sequence of length "lenSequence" for a CSR query,
 * 		  starting from the set of COIs passed to the application.
 *
 */
void UndirectedRN::setupQuerySequence(const uint32_t& lenSequence)
{
	// Case where we generate a sequence in the order specified by the user when it passes the COIs
	// via command line.
	if(lenSequence == 0)
	{
		this->startPosTypes = this->startPosTypesBackup;
		this->endPosTypes = this->endPosTypesBackup;
		this->sortedVecCostPOIs = this->sortedVecCostPOIsBackup;
		this->vecCostPOIs = this->vecCostPOIsBackup;
	}
	else
	{
		std::cout << "INIT => Setup COI sequence: ";
		this->startPosTypes.clear();
		this->endPosTypes.clear();
		this->sortedVecCostPOIs.clear();
		this->vecCostPOIs.clear();

		uint32_t lastIdx = rand() % this->sortedVecCostPOIsBackup.size(); std::cout << lastIdx << " ";
		this->sortedVecCostPOIs.push_back(this->sortedVecCostPOIsBackup[lastIdx]);
		this->startPosTypes.push_back(this->startPosTypesBackup[lastIdx]);
		this->endPosTypes.push_back(this->endPosTypesBackup[lastIdx]);
		this->vecCostPOIs.push_back(this->vecCostPOIsBackup[lastIdx]);
		for(uint32_t i = 1; i < lenSequence; i++)
		{
			uint32_t currIdx = rand() % this->sortedVecCostPOIsBackup.size();
			if(currIdx == lastIdx) currIdx = (currIdx + 1) % this->sortedVecCostPOIsBackup.size();
			lastIdx = currIdx;

			std::cout << lastIdx << " ";
			this->sortedVecCostPOIs.push_back(this->sortedVecCostPOIsBackup[lastIdx]);
			this->startPosTypes.push_back(this->startPosTypesBackup[lastIdx]);
			this->endPosTypes.push_back(this->endPosTypesBackup[lastIdx]);
			this->vecCostPOIs.push_back(this->vecCostPOIsBackup[lastIdx]);
		}
		std::cout << "\n";
	}
}

/**
 * @brief This method employs Dijkstra to compute the shortest path between a starting point "s"
 * 		  and a destination "d".
 *
 * @param s ID of the source.
 * @param d ID of the destination.
 *
 * @return Length of the shortest path.
 */
double UndirectedRN::SPDijkstra(const uint32_t& s, const uint32_t& d)
{
	// Instantiate the shortest path solver.
//	std::vector<Graph::Edge> tree_vector;
	lemon::Dijkstra<Graph, Graph::EdgeMap<double>> solver(this->g, this->mapEdges);

	// Find the shortest path between "s" and "d"
	solver.run(this->vecNodes[s], this->vecNodes[d]);

	// Retrieve the path computed by the solver.
	// lemon::Path<Graph> path = solver.path(this->vecNodes[d]);
	return(solver.dist(this->vecNodes[d]));


//	lemon::Path<Graph>::ArcIt it(path);
//	// NOTE: The iterator starts from the source and terminates at the destination.
//	std::cout << "Src: " << g.id(this->vecNodes[s]) << " -- Dest: " << g.id(this->vecNodes[d]) << "\n";
//	Graph::Arc arc;
//	double prova = 0;
//	while(g.target(it) != this->vecNodes[d])
//	{
//		arc = it;
//		prova += this->mapEdges[arc];
//		std::cout << "Edge " << g.id(it) << " (length " << this->mapEdges[arc] << ") goes from node "
//			<< g.id(g.source(it)) << " to node " << g.id(g.target(it)) << std::endl;
//		++it;
//	}
//	arc = it;
//	std::cout << "Edge " << g.id(it) << " goes from node "
//				<< g.id(g.source(it)) << " to node " << g.id(g.target(it)) << std::endl;
//	prova += this->mapEdges[arc];
//	std::cout << "Prova " << prova << "\n";
//
//
//	std::cout << "Euclidean distance points: " << this->computeDistanceVertices(this->vecNodes[s], this->vecNodes[d]) << " meters\n";
//	std::cout << "Cost of the path: " << solver.dist(this->vecNodes[d]) << " meters\n";
}

double UndirectedRN::SPAStar(const uint32_t& s, const uint32_t& d)
{
	typedef std::pair<double, uint32_t> ElQueue;
	std::priority_queue<ElQueue, std::vector<ElQueue>, CmpPQueue<ElQueue>> queueVertices;
	std::vector<bool> visitedVertices(this->vecNodes.size(), false);

	for(uint32_t i = 0; i < this->getNumVertices(); i++)
		queueVertices.push(ElQueue(0, i));


	// queueVertices.push(ElQueue(0 + this->computeDistancePointsIDs(s,d), s));
	queueVertices.push(ElQueue(0, s));
	ElQueue el;
	while(!queueVertices.empty())
	{
		// Pop the top element in the queue.
		el = queueVertices.top(); queueVertices.pop();
		// std::cout << "Popped vertex: " << el.second << "," << el.first <<"\n";

		// If this vertex was already marked as visited, discard it and pop the next
		// element in the queue.
		if(visitedVertices[el.second] == true) continue;

		// Mark this vertex as visited (we just found the shortest possible path touching it).
		visitedVertices[el.second] = true;

		// Check if the popped element is the destination.
		if(el.second == d) break;

		// 2 - Otherwise, find out the neighbors of this vertex and add them to the queue.
		const Graph::Node n = this->g.nodeFromId(el.second);
		// el.first -= this->computeDistancePointsIDs(el.second,d);
		for (Graph::IncEdgeIt e(this->g, n); e != lemon::INVALID; ++e)
		{
			// Check if this vertex was already visited.
			const uint32_t idCandidate = this->g.id(this->g.oppositeNode(n, static_cast<Graph::Edge>(e)));
			if(!visitedVertices[idCandidate])
			{
				// Add the vertex to the queue.
				queueVertices.push(std::pair<double, uint32_t>
			//					   (el.first + this->mapEdges[static_cast<Graph::Edge>(e)] + this->computeDistancePointsIDs(idCandidate,d), // Update weight path.
				   (el.first + this->mapEdges[static_cast<Graph::Edge>(e)], // Update weight path.
									idCandidate)); // ID neighbor.
			}
		}
	}

	return(el.first);
}

void UndirectedRN::precomputeEDPOIs(const std::string& outputFile)
{
	tbb::concurrent_unordered_map<uint32_t, tbb::concurrent_unordered_map<uint32_t, double>> tableEDs;

	// Calculate the Euclidean distances
	for(uint32_t category = 0; category < this->sortedVecCostPOIsBackup.size() - 1; category++)
	{
		for(uint32_t category2 = category + 1; category2 < this->sortedVecCostPOIsBackup.size(); category2++)
		{
			for(uint32_t p1 = 0; p1 < this->sortedVecCostPOIsBackup[category].size(); p1++)
			{
				const uint32_t source = this->sortedVecCostPOIsBackup[category][p1].second;
				for(uint32_t p2 = 0; p2 < this->sortedVecCostPOIsBackup[category2].size(); p2++)
				{
					const uint32_t target = this->sortedVecCostPOIsBackup[category2][p2].second;
					const double Edist = computeDistancePointsIDs(source, target);
					tableEDs[source][target] = Edist;
				}
			}
		}
	}

	// Write out the content of the hash table in a binary file.
	std::ofstream out(outputFile);
	for(const auto& p1 : tableEDs)
		for(const auto& p2 : p1.second)
			out << p1.first << " " << p2.first << " " << p2.second << "\n";
	out.close();

	// Write the number of POIs for each category (may be used to double check if the table can be used with
	// the COIs loaded in memory).
}

void UndirectedRN::precomputeTDPOIs(const std::string& outputFile)
{
	tbb::concurrent_unordered_map<uint32_t, tbb::concurrent_unordered_map<uint32_t, double>> tableTDs;

	// Calculate the travel distances.
	std::cout << "PRECOMPUTE => Using " << omp_get_max_threads() << " logical cores.\n";
	auto start_scoring = std::chrono::high_resolution_clock::now();
	uint32_t computedTDs = 0;
	for(uint32_t category = 0; category < this->sortedVecCostPOIsBackup.size() - 1; category++)
	{
		for(uint32_t category2 = category + 1; category2 < this->sortedVecCostPOIsBackup.size(); category2++)
		{
			for(uint32_t p1 = 0; p1 < this->sortedVecCostPOIsBackup[category].size(); p1++)
			{
				const uint32_t source = this->sortedVecCostPOIsBackup[category][p1].second;

				#pragma omp parallel for num_threads(omp_get_max_threads())
				for(uint32_t p2 = 0; p2 < this->sortedVecCostPOIsBackup[category2].size(); p2++)
				{
					const uint32_t target = this->sortedVecCostPOIsBackup[category2][p2].second;
					const double dist = this->SPDijkstra(source, target);
					tableTDs[source][target] = dist;
				}
				computedTDs += this->sortedVecCostPOIsBackup[category2].size();
				std::cout << "Computed SPs: " << computedTDs << "\n";
			}
		}
	}
	auto end_scoring = std::chrono::high_resolution_clock::now();
	auto total_TD_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_scoring - start_scoring).count();
	std::cout << "Time used to compute TDs: " << total_TD_time << " sec.\n";


	// Write out the content of the hash table in a binary file.
	std::ofstream out(outputFile);
	for(const auto& p1 : tableTDs)
		for(const auto& p2 : p1.second)
			out << p1.first << " " << p2.first << " " << p2.second << "\n";
	out.close();

	// Write the number of POIs for each category (may be used to double check if the table can be used with
	// the COIs loaded in memory).
}

void UndirectedRN::precomputeTDPOIsImproved(const std::string& outputFile)
{
	tbb::concurrent_unordered_map<uint32_t, tbb::concurrent_unordered_map<uint32_t, double>> tableTDs;

	// TODO: test the correctness.
	// Calculate the travel distances.
	std::cout << "PRECOMPUTE => Using " << omp_get_max_threads() << " logical cores.\n";
	auto start_scoring = std::chrono::high_resolution_clock::now();
	for(uint32_t category = 0; category < this->sortedVecCostPOIsBackup.size() - 1; category++)
	{
		std::unordered_set<uint32_t> setTarget;
		for(uint32_t category2 = category + 1; category2 < this->sortedVecCostPOIsBackup.size(); category2++)
		{
			for(const auto& p : this->sortedVecCostPOIsBackup[category2])
				setTarget.insert(p.second);
		}

		std::cout << "Computing distances from POIs of category " << category << "\n";
		#pragma omp parallel for num_threads(omp_get_max_threads())
		for(uint32_t p1 = 0; p1 < this->sortedVecCostPOIsBackup[category].size(); p1++)
		{
			const uint32_t& source = this->sortedVecCostPOIsBackup[category][p1].second;
			this->findNNDijkstra(source, setTarget, tableTDs);
		}
	}
	auto end_scoring = std::chrono::high_resolution_clock::now();
	auto total_TD_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_scoring - start_scoring).count();
	std::cout << "Time took to compute TDs: " << total_TD_time << " sec.\n";


	// Write out the content of the hash table in a binary file.
	std::ofstream out(outputFile);
	for(const auto& p1 : tableTDs)
		for(const auto& p2 : p1.second)
			out << p1.first << " " << p2.first << " " << p2.second << "\n";
	out.close();

	// Write the number of POIs for each category (may be used to double check if the table can be used with
	// the COIs loaded in memory).
}

void UndirectedRN::findNNDijkstra(const uint32_t& src,
								  const std::unordered_set<uint32_t>& targetSet,
								  tbb::concurrent_unordered_map<uint32_t, tbb::concurrent_unordered_map<uint32_t, double>>& tableTDs)
{
	typedef std::pair<double, uint32_t> ElQueue;
	std::priority_queue<ElQueue, std::vector<ElQueue>, CmpPQueue<ElQueue>> queueVertices;
	std::vector<bool> visitedVertices(this->vecNodes.size(), false);

	queueVertices.push(ElQueue(0, src));
	uint32_t cnt = 0;
	ElQueue el;
	while(!queueVertices.empty())
	{
		// Pop the top element in the queue.
		el = queueVertices.top(); queueVertices.pop();
		// std::cout << "Popped vertex: " << el.second << "," << el.first <<"\n";

		// If this vertex was already marked as visited, discard it and pop the next
		// element in the queue.
		if(visitedVertices[el.second] == true) continue;

		// Mark this vertex as visited (we just found the shortest possible path touching it).
		visitedVertices[el.second] = true;

		// 1 - Check if the popped element is associated with a POI: if so, return the ID of the POI.
		if(targetSet.find(el.second) != targetSet.end())
		{
			tableTDs[src][el.second] = el.first;
			if(++cnt == targetSet.size()) return;
		}

		// 2 - Otherwise, find out the neighbors of this vertex and add them to the queue.
		const Graph::Node n = this->g.nodeFromId(el.second);
		for (Graph::IncEdgeIt e(this->g, n); e != lemon::INVALID; ++e)
		{
			// Check if the currently considered neighbor was already visited before.
			const uint32_t idCandidate = this->g.id(this->g.oppositeNode(n, static_cast<Graph::Edge>(e)));
			if(!visitedVertices[idCandidate])
			{
				// If not, add the vertex to the queue.
				queueVertices.push(std::pair<double, uint32_t>
								   (el.first + this->mapEdges[static_cast<Graph::Edge>(e)], // Update weight path.
									idCandidate)); // ID neighbor.

	//		    std::cout << "Edge " << this->g.id(static_cast<Graph::Edge>(e)) << "\n";
//				std::cout << "Vertex queued: " << this->g.id(this->g.oppositeNode(this->g.nodeFromId(el.second), static_cast<Graph::Edge>(e))) <<
//							 " - distance: " << el.first + this->mapEdges[static_cast<Graph::Edge>(e)] << "\n";
			}
		}
	}
}

void UndirectedRN::loadPrecomputedTDPOIs(const std::string& nameFile)
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
		if(this->tableTDs.find(source) == this->tableTDs.end())
			this->tableTDs[source] = std::vector<std::pair<double, uint32_t>>();
		this->tableTDs[source].push_back(std::make_pair(dist, target));

		if(this->tableTDs.find(target) == this->tableTDs.end())
			this->tableTDs[target] = std::vector<std::pair<double, uint32_t>>();
		this->tableTDs[target].push_back(std::make_pair(dist, source));
	}
	infile.close();

	// Sort the POIs associated with each POI in increasing order of travel distance.
	for(auto& poi : this->tableTDs)
		std::sort(poi.second.begin(), poi.second.end());
}

void UndirectedRN::initPrecomputedNNsPOIs()
{
	// For each POI in the table, create another table where for each COI we have the associated POIs
	// in increasing order of travel distance w.r.t. the considered poi.
	const uint32_t maxTypes = this->startPosTypes.size() - 1;
	for(uint32_t type = 0; type < maxTypes; type++)
	{
		for(uint32_t idPOI = this->startPosTypes[type]; idPOI <= this->endPosTypes[type]; idPOI++)
		{
			// Create the entry in the table, if needed.
			if(this->tableTypeTDs.find(idPOI) == this->tableTypeTDs.end())
				this->tableTypeTDs[idPOI] = std::vector<std::vector<std::pair<double, uint32_t>>>(sortedVecCostPOIs.size());

			// Scan the POIs associated with idPOI, and insert the ones belonging to type "type+1".
			for(const auto& el : this->tableTDs[idPOI])
			{
				// Add the POI to the vector.
				if(el.second >= this->startPosTypes[type+1] && el.second <= this->endPosTypes[type+1])
				{
					this->tableTypeTDs[idPOI][type+1].push_back(std::make_pair(el.first, el.second));
				}
			}
		}
	}
}

/**
 * @brief This method implements the Incremental Network Expansion.
 *
 * @param idSource ID of the vertex representing the starting point.
 * @param typePOI Category of the POi we are targeting.
 *
 * @return ID of NN vertex belonging to the category of interest.
 */
std::pair<int32_t, double> UndirectedRN::INELookup(const uint32_t& idSource, const uint32_t& typePOI, const uint32_t& k)
{
	// 1 - Check if the k-NN is available in the lookup table.
	if(this->tableTypeTDs.find(idSource) != this->tableTypeTDs.end())
	{
		if(k < this->tableTypeTDs[idSource][typePOI].size())
		{
			const auto& kNN = this->tableTypeTDs[idSource][typePOI][k];
			return(std::make_pair(kNN.second, kNN.first));
		}
	}
	else
	{
		this->tableTypeTDs[idSource] =
				std::vector<std::vector<std::pair<double, uint32_t>>>(this->startPosTypes.size());
	}


	// 2 - If the k-NN is not available in the lookup table, then proceed with INE.
	typedef std::pair<double, uint32_t> ElQueue;
	std::priority_queue<ElQueue, std::vector<ElQueue>, CmpPQueue<ElQueue>> queueVertices;
	std::vector<bool> visitedVertices(this->vecNodes.size(), false);


//	std::cout << "Finding NN POI of category " << typePOI << " for vertex " << idSource << "\n";
//	std::cout << "Range: [" << startPosTypes[typePOI] << "," << endPosTypes[typePOI] << "]\n";

	queueVertices.push(ElQueue(0, idSource));
	ElQueue el;
	uint32_t cnt = 0;
	while(!queueVertices.empty())
	{
		// Pop the top element in the queue.
		el = queueVertices.top(); queueVertices.pop();
		// std::cout << "Popped vertex: " << el.second << "," << el.first <<"\n";

		// If this vertex was already marked as visited, discard it and pop the next
		// element in the queue.
		if(visitedVertices[el.second] == true) continue;

		// Mark this vertex as visited (we just found the shortest possible path touching it).
		visitedVertices[el.second] = true;

		// 1 - Check if the popped element is associated with a POI: if so, return the ID of the POI.
		if(el.second >= startPosTypes[typePOI] && el.second <= endPosTypes[typePOI])
		{
			// Check if we found out again a POI that we have to ignore.
			if(cnt++ == k) break;
		}

		// 2 - Otherwise, find out the neighbors of this vertex and add them to the queue.
		const Graph::Node n = this->g.nodeFromId(el.second);
		for (Graph::IncEdgeIt e(this->g, n); e != lemon::INVALID; ++e)
		{
			// Check if this vertex was already visited.
			const uint32_t idCandidate = this->g.id(this->g.oppositeNode(n, static_cast<Graph::Edge>(e)));
			if(!visitedVertices[idCandidate])
			{
				// Add the vertex to the queue.
				queueVertices.push(std::pair<double, uint32_t>
								   (el.first + this->mapEdges[static_cast<Graph::Edge>(e)], // Update weight path.
									idCandidate)); // ID neighbor.

	//		    std::cout << "Edge " << this->g.id(static_cast<Graph::Edge>(e)) << "\n";
//				std::cout << "Vertex queued: " << this->g.id(this->g.oppositeNode(this->g.nodeFromId(el.second), static_cast<Graph::Edge>(e))) <<
//							 " - distance: " << el.first + this->mapEdges[static_cast<Graph::Edge>(e)] << "\n";
			}
		}
	}


	this->tableTypeTDs[idSource][typePOI].push_back(el);
	return(std::pair<int32_t, double>(el.second, el.first));
}

/**
 * @brief This method implements the Incremental Network Expansion.
 *
 * @param idSource ID of the vertex representing the starting point.
 * @param typePOI Category of the POi we are targeting.
 *
 * @return ID of NN vertex belonging to the category of interest.
 */
void UndirectedRN::initkNNsSource(const uint32_t& idSource)
{
	// Initialize entry source in hash table (ONLY IF NEEDED, the source may be itself a POI so
	// a new allocation would destroy the information about kNNs in other COIs!).
	if(this->tableTypeTDs.find(idSource) == this->tableTypeTDs.end())
		this->tableTypeTDs[idSource] =
				std::vector<std::vector<std::pair<double, uint32_t>>>(this->startPosTypes.size());


	typedef std::pair<double, uint32_t> ElQueue;
	std::priority_queue<ElQueue, std::vector<ElQueue>, CmpPQueue<ElQueue>> queueVertices;
	std::vector<bool> visitedVertices(this->vecNodes.size(), false);


//	std::cout << "Finding NN POI of category " << typePOI << " for vertex " << idSource << "\n";
//	std::cout << "Range: [" << startPosTypes[typePOI] << "," << endPosTypes[typePOI] << "]\n";

	queueVertices.push(ElQueue(0, idSource));
	ElQueue el;
	uint32_t k = 0;
	const uint32_t numPOIs = (this->endPosTypes[0] - this->startPosTypes[0]) + 1;
	while(!queueVertices.empty())
	{
		// Pop the top element in the queue.
		el = queueVertices.top(); queueVertices.pop();
		// std::cout << "Popped vertex: " << el.second << "," << el.first <<"\n";

		// If this vertex was already marked as visited, discard it and pop the next
		// element in the queue.
		if(visitedVertices[el.second] == true) continue;

		// Mark this vertex as visited (we just found the shortest possible path touching it).
		visitedVertices[el.second] = true;

		// 1 - Check if the popped element is associated with a POI of the COI we are considering.
		if(el.second >= this->startPosTypes[0] && el.second <= this->endPosTypes[0])
		{
			// Check if we found out again a POI that we have to ignore.
			this->tableTypeTDs[idSource][0].push_back(el);
			if(++k == numPOIs) break;
		}

		// 2 - Otherwise, find out the neighbors of this vertex and add them to the queue.
		const Graph::Node n = this->g.nodeFromId(el.second);
		for (Graph::IncEdgeIt e(this->g, n); e != lemon::INVALID; ++e)
		{
			// Check if this vertex was already visited.
			const uint32_t idCandidate = this->g.id(this->g.oppositeNode(n, static_cast<Graph::Edge>(e)));
			if(!visitedVertices[idCandidate])
			{
				// Add the vertex to the queue.
				queueVertices.push(std::pair<double, uint32_t>
								   (el.first + this->mapEdges[static_cast<Graph::Edge>(e)], // Update weight path.
									idCandidate)); // ID neighbor.

	//		    std::cout << "Edge " << this->g.id(static_cast<Graph::Edge>(e)) << "\n";
//				std::cout << "Vertex queued: " << this->g.id(this->g.oppositeNode(this->g.nodeFromId(el.second), static_cast<Graph::Edge>(e))) <<
//							 " - distance: " << el.first + this->mapEdges[static_cast<Graph::Edge>(e)] << "\n";
			}
		}
	}
}

/**
 * @brief This method implements the Progressive Network Expansion to compute the optimal sequenced route.
 *
 * @param idSource ID of the starting vertex.
 */
void UndirectedRN::PNEBaseline(const uint32_t& idSource, const double& minCost)
{
	// Some typedefs and declarations.
	typedef struct
	{   std::vector<uint32_t> elSeq; std::vector<uint32_t> elK; std::vector<double> elW;
	} Sequence;
	typedef std::pair<double, Sequence> ElQueue;
	typedef SimpleSkyline::ElementSkyline ElementSkyline;
	std::priority_queue<ElQueue, std::vector<ElQueue>, CmpPQueue<ElQueue>> queuePaths;


	// Explore the search space.
	const uint32_t numTypes = this->sortedVecCostPOIs.size();
	queuePaths.push(ElQueue(0, {std::vector<uint32_t>({idSource}), std::vector<uint32_t>({0}), std::vector<double>({0})}));
	ElQueue el;
	double POICostUB = std::numeric_limits<double>::max();
	SimpleSkyline* LS = 0;
	while(!queuePaths.empty())
	{
		// Pop the top element of the queue.
		el = queuePaths.top(); queuePaths.pop();

		double POIcost = 0;
		for(uint32_t i = 1; i < el.second.elSeq.size(); i++)
			POIcost += this->getCostPOI(i-1, el.second.elSeq[i]);


		// DEBUG.
//		std::cout << "Popped sequence ";
//		for(auto elem : el.second.elSeq)
//			std::cout << elem << " ";
//		std::cout << "\nPopped sequence cost ";
//		for(uint32_t i = 1; i < el.second.elSeq.size(); i++)
//			std::cout << this->getCostPOI(i-1, el.second.elSeq[i]) << " ";
//		std::cout << "\n";
//		std::cout << "(W " << el.first << ", C "<< POIcost << ")\n";
//		std::cout << "TCUB: " << POICostUB << "\n";
//		std::cout << "TDUB: " << ubTD << "\n";

		// If we pass the travel distance upper bound, terminate the search.
//		if(el.first > ubTD)
//		{
////			std::cout << "Early terminate!\n";
//			break;
//		}

		// 1 - Check if we have found a full path: if so, attempt to insert the sequence into the skyline.
		if(el.second.elSeq.size() == numTypes + 1)
		{
//			std::cout << "Found full route\n";
			// Case where we need to initialize the linear skyline.
			if(LS == 0)
			{
				LS = new SimpleSkyline(ElementSkyline(el.first, POIcost), el.second.elSeq);
				POICostUB = POIcost;

				// DEBUG.
//				std::cout << "Init skyline with R^{TD}! ";
//				for(auto elem : el.second.elSeq)
//					std::cout << elem << " ";
//				std::cout << "(W " << el.first << ", C "<< POIcost << ")\n";
			}

			// Check if we can insert the route into the skyline.
			if(POIcost <= POICostUB)
			{
				// Attempt to insert the route into the skyline.
				if(LS->insertElement(ElementSkyline(el.first, POIcost), el.second.elSeq))
				{
					POICostUB = POIcost;

					// DEBUG.
//					std::cout << "Found non dominated route! ";
//					for(auto elem : el.second.elSeq)
//						std::cout << elem << " ";
//					std::cout << "(W " << el.first << ", C "<< POIcost << ")\n";
				}

				// Terminate the algorithm when we find the first route
				// with the minimum possible travel cost.
				if(POIcost == minCost)
				{
					// ubTD = el.first;
					// std::cout << "EARLY\n";
					break;
				}
			}
		}
		// 2 - Otherwise, if the sequence is below the POI cost upper bound,
		//     extend it with the NN of the last COI in the sequence.
		else
		{
			if(POIcost <= POICostUB)
			{
				const auto NN = this->INELookup(el.second.elSeq.back(),		// Source node.
						el.second.elSeq.size() - 1,   					    // Target category.
						el.second.elK.back());							    // Find k-th NN.


				// 3 - Generate the new partial sequence;
				Sequence tmpSeq = el.second;
				tmpSeq.elSeq.push_back(NN.first); tmpSeq.elW.push_back(tmpSeq.elW.back() + NN.second); tmpSeq.elK.push_back(0);
				queuePaths.push(ElQueue(el.first + NN.second, tmpSeq));

				// DEBUG.
//							std::cout << "Computing NN\n";
				//			std::cout << "Generated sequence: ";
				//			for(auto elem : tmpSeq.elSeq)
				//				std::cout << elem << " ";
				//			std::cout << "(W " << el.first + NN.second << ")\n";
			}
//			else
//				std::cout << "Not extending popped sequence\n";
		}


		// 4 - Look for the next NN for the predecessor of el. Then, add it to the queue.
		if(el.second.elSeq.size() > 1)
		{
			Sequence& seq = el.second;

			seq.elSeq.pop_back(); seq.elK.pop_back(); seq.elW.pop_back(); // Remove the last category from the sequence.
			seq.elK.back()++; 		   // Update the k of the next NN to find.
			el.first = seq.elW.back(); // Update the partial weight in el accordingly.

			double POIcostPred = 0;
			for(uint32_t i = 1; i < seq.elSeq.size(); i++)
				POIcostPred += this->getCostPOI(i-1, seq.elSeq[i]);
//			std::cout << "Attempt to generate prec. seq.: " << POIcostPred << "\n";

			// Find the k-th NN (if possible).
//			const auto kappa = seq.elK.back();
//			std::cout << "TEST: " << seq.elK.size() - 1 << "," << this->sortedVecCostPOIs[seq.elK.size() - 1].size() << "\n";
			if(POIcostPred <= POICostUB && seq.elK.back() < this->sortedVecCostPOIs[seq.elK.size() - 1].size())
			{
//				std::cout << "Computing NN2\n";
				const auto NNPred = this->INELookup(seq.elSeq.back(),		// Source node.
											  seq.elSeq.size() - 1, // Target category.
											  seq.elK.back());		// Find k-th NN.

				// Update el and push it in the queue.
				seq.elSeq.push_back(NNPred.first); seq.elW.push_back(seq.elW.back() + NNPred.second); seq.elK.push_back(0);
				el.first = seq.elW.back();
				queuePaths.push(el);

				// DEBUG.
//				std::cout << "Generated prec. sequence (k=" << kappa << "): ";
//				for(auto elem : el.second.elSeq)
//					std::cout << elem << " ";
//				std::cout << "(W " << el.first << ")\n";
			}
//			else
//				std::cout << "Reached END!!\n";
		}
	}


	std::cout << "Number of elements in the skyline: " << LS->getLS().size() << "\n";
	for(auto els : LS->getLS())
		std::cout << els.first << "," << els.second << "\n";
	for(auto els : LS->getLSS())
	{
		for(auto elp : els)
			std::cout << elp << " ";
		std::cout << "\n";
	}
}

/**
 * @brief This method substitutes the Incremental Network Expansion when the distance between POIs
 * 		  are precomputed.
 *
 * @param idSource ID of the vertex representing the starting point.
 * @param typePOI Category of the POi we are targeting.
 *
 * @return ID of NN vertex belonging to the category of interest.
 */
std::pair<int32_t, double> UndirectedRN::INEPrecomp(const uint32_t& idSource, const uint32_t& typePOI, const uint32_t& k)
{
	const auto& dest = this->tableTypeTDs[idSource][typePOI][k];
	return(std::pair<int32_t, double>(dest.second, dest.first));
}

/**
 * @brief This method implements the Progressive Network Expansion to compute the optimal sequenced route.
 *
 * @param idSource ID of the starting vertex.
 */
void UndirectedRN::PNEBaselinePrecomp(const uint32_t& idSource, const double& minCost)
{
	static uint32_t numElsSkyline = 0;
	static uint32_t numRounds = 0;


	// Some typedefs and declarations.
	typedef struct
	{   std::vector<uint32_t> elSeq; std::vector<uint32_t> elK; std::vector<double> elW;
	} Sequence;
	typedef std::pair<double, Sequence> ElQueue;
	typedef SimpleSkyline::ElementSkyline ElementSkyline;
	std::priority_queue<ElQueue, std::vector<ElQueue>, CmpPQueue<ElQueue>> queuePaths;


	// Explore the search space.
	const uint32_t numTypes = this->sortedVecCostPOIs.size();
	queuePaths.push(ElQueue(0, {std::vector<uint32_t>({idSource}), std::vector<uint32_t>({0}), std::vector<double>({0})}));
	ElQueue el;
	double POICostUB = std::numeric_limits<double>::max();
	SimpleSkyline* LS = 0;
	while(!queuePaths.empty())
	{
		// Pop the top element of the queue.
		el = queuePaths.top(); queuePaths.pop();

		double POIcost = 0;
		for(uint32_t i = 1; i < el.second.elSeq.size(); i++)
			POIcost += this->getCostPOI(i-1, el.second.elSeq[i]);


		// DEBUG.
//		std::cout << "Popped sequence ";
//		for(auto elem : el.second.elSeq)
//			std::cout << elem << " ";
//		std::cout << "\nPopped sequence cost ";
//		for(uint32_t i = 1; i < el.second.elSeq.size(); i++)
//			std::cout << this->getCostPOI(i-1, el.second.elSeq[i]) << " ";
//		std::cout << "\n";
//		std::cout << "(W " << el.first << ", C "<< POIcost << ")\n";
//		std::cout << "TCUB: " << POICostUB << "\n";
//		std::cout << "TDUB: " << ubTD << "\n";

		// If we pass the travel distance upper bound, terminate the search.
//		if(el.first > ubTD)
//		{
////			std::cout << "Early terminate!\n";
//			break;
//		}

		// 1 - Check if we have found a full path: if so, attempt to insert the sequence into the skyline.
		if(el.second.elSeq.size() == numTypes + 1)
		{
//			std::cout << "Found full route\n";
			// Case where we need to initialize the linear skyline.
			if(LS == 0)
			{
				LS = new SimpleSkyline(ElementSkyline(el.first, POIcost), el.second.elSeq);
				POICostUB = POIcost;

				// DEBUG.
//				std::cout << "Init skyline with R^{TD}! ";
//				for(auto elem : el.second.elSeq)
//					std::cout << elem << " ";
//				std::cout << "(W " << el.first << ", C "<< POIcost << ")\n";
			}

			// Check if we can insert the route into the skyline.
			if(POIcost <= POICostUB)
			{
				// Attempt to insert the route into the skyline.
				if(LS->insertElement(ElementSkyline(el.first, POIcost), el.second.elSeq))
				{
					POICostUB = POIcost;

					// DEBUG.
//					std::cout << "Found non dominated route! ";
//					for(auto elem : el.second.elSeq)
//						std::cout << elem << " ";
//					std::cout << "(W " << el.first << ", C "<< POIcost << ")\n";
				}
			}

			// Terminate the algorithm when we find the first route
			// with the minimum possible travel cost.
			if(POIcost == minCost)
			{
				// ubTD = el.first;
				// std::cout << "EARLY\n";
				break;
			}
		}
		// 2 - Otherwise, if the sequence is below the POI cost upper bound,
		//     extend it with the NN of the last COI in the sequence.
		else
		{
			if(POIcost <= POICostUB)
			{
				const auto NN = (el.second.elSeq.size() - 1) == 0 ?
									  this->INELookup(el.second.elSeq.back(),	// Source node.
									  el.second.elSeq.size() - 1,   	// Target category.
									  el.second.elK.back()) :			// Find k-th NN.
									  this->INEPrecomp(el.second.elSeq.back(),		// Source node.
									  el.second.elSeq.size() - 1,   			    // Target category.
									  el.second.elK.back());					    // Find k-th NN.


				// 3 - Generate the new partial sequence;
				Sequence tmpSeq = el.second;
				tmpSeq.elSeq.push_back(NN.first); tmpSeq.elW.push_back(tmpSeq.elW.back() + NN.second); tmpSeq.elK.push_back(0);
				queuePaths.push(ElQueue(el.first + NN.second, tmpSeq));

				// DEBUG.
//							std::cout << "Computing NN\n";
				//			std::cout << "Generated sequence: ";
				//			for(auto elem : tmpSeq.elSeq)
				//				std::cout << elem << " ";
				//			std::cout << "(W " << el.first + NN.second << ")\n";
			}
//			else
//				std::cout << "Not extending popped sequence\n";
		}


		// 4 - Look for the next NN for the predecessor of el. Then, add it to the queue.
		if(el.second.elSeq.size() > 1)
		{
			Sequence& seq = el.second;

			seq.elSeq.pop_back(); seq.elK.pop_back(); seq.elW.pop_back(); // Remove the last category from the sequence.
			seq.elK.back()++; 		   // Update the k of the next NN to find.
			el.first = seq.elW.back(); // Update the partial weight in el accordingly.

			double POIcostPred = 0;
			for(uint32_t i = 1; i < seq.elSeq.size(); i++)
				POIcostPred += this->getCostPOI(i-1, seq.elSeq[i]);
//			std::cout << "Attempt to generate prec. seq.: " << POIcostPred << "\n";

			// Find the k-th NN (if possible).
//			const auto kappa = seq.elK.back();
//			std::cout << "TEST: " << seq.elK.size() - 1 << "," << this->sortedVecCostPOIs[seq.elK.size() - 1].size() << "\n";
			if(POIcostPred <= POICostUB && seq.elK.back() < this->sortedVecCostPOIs[seq.elK.size() - 1].size())
			{
//				std::cout << "Computing NN2\n";
				const auto NNPred = (seq.elSeq.size() - 1 == 0) ?
										this->INELookup(seq.elSeq.back(),	// Source node.
										seq.elSeq.size() - 1,  	    	    // Target category.
										seq.elK.back()) :					// Find k-th NN.
										this->INEPrecomp(seq.elSeq.back(),	// Source node.
										seq.elSeq.size() - 1,  	    		// Target category.
										seq.elK.back());					// Find k-th NN.


				// Update el and push it in the queue.
				seq.elSeq.push_back(NNPred.first); seq.elW.push_back(seq.elW.back() + NNPred.second); seq.elK.push_back(0);
				el.first = seq.elW.back();
				queuePaths.push(el);

				// DEBUG.
//				std::cout << "Generated prec. sequence (k=" << kappa << "): ";
//				for(auto elem : el.second.elSeq)
//					std::cout << elem << " ";
//				std::cout << "(W " << el.first << ")\n";
			}
//			else
//				std::cout << "Reached END!!\n";
		}
	}


	std::cout << "Number of elements in the skyline: " << LS->getLS().size() << "\n";
	for(auto els : LS->getLS())
		std::cout << els.first << "," << els.second << "\n";
	for(auto els : LS->getLSS())
	{
		for(auto elp : els)
			std::cout << elp << " ";
		std::cout << "\n";
	}


	numElsSkyline += LS->getLS().size(); numRounds++;
	std::cout << "Average number of elements in the skyline: " << ((double) numElsSkyline) / numRounds << "\n";
}

/**
 * @brief This method implements the Progressive Network Expansion to compute the optimal sequenced route.
 *
 * @param idSource ID of the starting vertex.
 */
void UndirectedRN::PNEBaselinePrecompConvSkyline(const uint32_t& idSource, const double& minCost)
{
	static uint32_t numElsSkyline = 0;
	static uint32_t numRounds = 0;


	// Some typedefs and declarations.
	typedef struct
	{   std::vector<uint32_t> elSeq; std::vector<uint32_t> elK; std::vector<double> elW;
	} Sequence;
	typedef std::pair<double, Sequence> ElQueue;
	typedef SimpleSkyline::ElementSkyline ElementSkyline;
	std::priority_queue<ElQueue, std::vector<ElQueue>, CmpPQueue<ElQueue>> queuePaths;


	// Explore the search space.
	const uint32_t numTypes = this->sortedVecCostPOIs.size();
	queuePaths.push(ElQueue(0, {std::vector<uint32_t>({idSource}), std::vector<uint32_t>({0}), std::vector<double>({0})}));
	ElQueue el;
	double POICostUB = std::numeric_limits<double>::max();
	SimpleSkyline* LS = 0;
	while(!queuePaths.empty())
	{
		// Pop the top element of the queue.
		el = queuePaths.top(); queuePaths.pop();

		double POIcost = 0;
		for(uint32_t i = 1; i < el.second.elSeq.size(); i++)
			POIcost += this->getCostPOI(i-1, el.second.elSeq[i]);


		// DEBUG.
//		std::cout << "Popped sequence ";
//		for(auto elem : el.second.elSeq)
//			std::cout << elem << " ";
//		std::cout << "\nPopped sequence cost ";
//		for(uint32_t i = 1; i < el.second.elSeq.size(); i++)
//			std::cout << this->getCostPOI(i-1, el.second.elSeq[i]) << " ";
//		std::cout << "\n";
//		std::cout << "(W " << el.first << ", C "<< POIcost << ")\n";
//		std::cout << "TCUB: " << POICostUB << "\n";
//		std::cout << "TDUB: " << ubTD << "\n";

		// If we pass the travel distance upper bound, terminate the search.
//		if(el.first > ubTD)
//		{
////			std::cout << "Early terminate!\n";
//			break;
//		}

		// 1 - Check if we have found a full path: if so, attempt to insert the sequence into the skyline.
		if(el.second.elSeq.size() == numTypes + 1)
		{
//			std::cout << "Found full route\n";
			// Case where we need to initialize the linear skyline.
			if(LS == 0)
			{
				LS = new SimpleConventionalSkyline(ElementSkyline(el.first, POIcost), el.second.elSeq);
				POICostUB = POIcost;

				// DEBUG.
//				std::cout << "Init skyline with R^{TD}! ";
//				for(auto elem : el.second.elSeq)
//					std::cout << elem << " ";
//				std::cout << "(W " << el.first << ", C "<< POIcost << ")\n";
			}

			// Check if we can insert the route into the skyline.
			if(POIcost <= POICostUB)
			{
				// Attempt to insert the route into the skyline.
				if(LS->insertElement(ElementSkyline(el.first, POIcost), el.second.elSeq))
				{
					POICostUB = POIcost;

					// DEBUG.
//					std::cout << "Found non dominated route! ";
//					for(auto elem : el.second.elSeq)
//						std::cout << elem << " ";
//					std::cout << "(W " << el.first << ", C "<< POIcost << ")\n";
				}
			}

			// Terminate the algorithm when we find the first route
			// with the minimum possible travel cost.
			if(POIcost == minCost)
			{
				// ubTD = el.first;
				// std::cout << "EARLY\n";
				break;
			}
		}
		// 2 - Otherwise, if the sequence is below the POI cost upper bound,
		//     extend it with the NN of the last COI in the sequence.
		else
		{
			if(POIcost <= POICostUB)
			{
				const auto NN = (el.second.elSeq.size() - 1) == 0 ?
									  this->INELookup(el.second.elSeq.back(),	// Source node.
									  el.second.elSeq.size() - 1,   	// Target category.
									  el.second.elK.back()) :			// Find k-th NN.
									  this->INEPrecomp(el.second.elSeq.back(),		// Source node.
									  el.second.elSeq.size() - 1,   			    // Target category.
									  el.second.elK.back());					    // Find k-th NN.


				// 3 - Generate the new partial sequence;
				Sequence tmpSeq = el.second;
				tmpSeq.elSeq.push_back(NN.first); tmpSeq.elW.push_back(tmpSeq.elW.back() + NN.second); tmpSeq.elK.push_back(0);
				queuePaths.push(ElQueue(el.first + NN.second, tmpSeq));

				// DEBUG.
//							std::cout << "Computing NN\n";
				//			std::cout << "Generated sequence: ";
				//			for(auto elem : tmpSeq.elSeq)
				//				std::cout << elem << " ";
				//			std::cout << "(W " << el.first + NN.second << ")\n";
			}
//			else
//				std::cout << "Not extending popped sequence\n";
		}


		// 4 - Look for the next NN for the predecessor of el. Then, add it to the queue.
		if(el.second.elSeq.size() > 1)
		{
			Sequence& seq = el.second;

			seq.elSeq.pop_back(); seq.elK.pop_back(); seq.elW.pop_back(); // Remove the last category from the sequence.
			seq.elK.back()++; 		   // Update the k of the next NN to find.
			el.first = seq.elW.back(); // Update the partial weight in el accordingly.

			double POIcostPred = 0;
			for(uint32_t i = 1; i < seq.elSeq.size(); i++)
				POIcostPred += this->getCostPOI(i-1, seq.elSeq[i]);
//			std::cout << "Attempt to generate prec. seq.: " << POIcostPred << "\n";

			// Find the k-th NN (if possible).
//			const auto kappa = seq.elK.back();
//			std::cout << "TEST: " << seq.elK.size() - 1 << "," << this->sortedVecCostPOIs[seq.elK.size() - 1].size() << "\n";
			if(POIcostPred <= POICostUB && seq.elK.back() < this->sortedVecCostPOIs[seq.elK.size() - 1].size())
			{
//				std::cout << "Computing NN2\n";
				const auto NNPred = (seq.elSeq.size() - 1 == 0) ?
										this->INELookup(seq.elSeq.back(),	// Source node.
										seq.elSeq.size() - 1,  	    	    // Target category.
										seq.elK.back()) :					// Find k-th NN.
										this->INEPrecomp(seq.elSeq.back(),	// Source node.
										seq.elSeq.size() - 1,  	    		// Target category.
										seq.elK.back());					// Find k-th NN.


				// Update el and push it in the queue.
				seq.elSeq.push_back(NNPred.first); seq.elW.push_back(seq.elW.back() + NNPred.second); seq.elK.push_back(0);
				el.first = seq.elW.back();
				queuePaths.push(el);

				// DEBUG.
//				std::cout << "Generated prec. sequence (k=" << kappa << "): ";
//				for(auto elem : el.second.elSeq)
//					std::cout << elem << " ";
//				std::cout << "(W " << el.first << ")\n";
			}
//			else
//				std::cout << "Reached END!!\n";
		}
	}


	std::cout << "Number of elements in the skyline: " << LS->getLS().size() << "\n";
	for(auto els : LS->getLS())
		std::cout << els.first << "," << els.second << "\n";
	for(auto els : LS->getLSS())
	{
		for(auto elp : els)
			std::cout << elp << " ";
		std::cout << "\n";
	}


	numElsSkyline += LS->getLS().size(); numRounds++;
	std::cout << "Average number of elements in the skyline: " << ((double) numElsSkyline) / numRounds << "\n";
}

/**
 * @brief This method implements the Progressive Network Expansion to compute the optimal sequenced route.
 *
 * @param idSource ID of the starting vertex.
 */
void UndirectedRN::PNEBaselinePrecompWithoutPruning(const uint32_t& idSource, const double& minCost)
{
	// Some typedefs and declarations.
	typedef struct
	{   std::vector<uint32_t> elSeq; std::vector<uint32_t> elK; std::vector<double> elW;
	} Sequence;
	typedef std::pair<double, Sequence> ElQueue;
	typedef SimpleSkyline::ElementSkyline ElementSkyline;
	std::priority_queue<ElQueue, std::vector<ElQueue>, CmpPQueue<ElQueue>> queuePaths;


	// Explore the search space.
	const uint32_t numTypes = this->sortedVecCostPOIs.size();
	queuePaths.push(ElQueue(0, {std::vector<uint32_t>({idSource}), std::vector<uint32_t>({0}), std::vector<double>({0})}));
	ElQueue el;
	double POICostUB = std::numeric_limits<double>::max();
	SimpleSkyline* LS = 0;
	while(!queuePaths.empty())
	{
		// Pop the top element of the queue.
		el = queuePaths.top(); queuePaths.pop();

		double POIcost = 0;
		for(uint32_t i = 1; i < el.second.elSeq.size(); i++)
			POIcost += this->getCostPOI(i-1, el.second.elSeq[i]);


		// DEBUG.
//		std::cout << "Popped sequence ";
//		for(auto elem : el.second.elSeq)
//			std::cout << elem << " ";
//		std::cout << "\nPopped sequence cost ";
//		for(uint32_t i = 1; i < el.second.elSeq.size(); i++)
//			std::cout << this->getCostPOI(i-1, el.second.elSeq[i]) << " ";
//		std::cout << "\n";
//		std::cout << "(W " << el.first << ", C "<< POIcost << ")\n";
//		std::cout << "TCUB: " << POICostUB << "\n";
//		std::cout << "TDUB: " << ubTD << "\n";


		// 1 - Check if we have found a full path: if so, attempt to insert the sequence into the skyline.
		if(el.second.elSeq.size() == numTypes + 1)
		{
//			std::cout << "Found full route\n";
			// Case where we need to initialize the linear skyline.
			if(LS == 0)
			{
				LS = new SimpleSkyline(ElementSkyline(el.first, POIcost), el.second.elSeq);
				POICostUB = POIcost;

				// DEBUG.
//				std::cout << "Init skyline with R^{TD}! ";
//				for(auto elem : el.second.elSeq)
//					std::cout << elem << " ";
//				std::cout << "(W " << el.first << ", C "<< POIcost << ")\n";
			}

			// Check if we can insert the route into the skyline.
			if(POIcost <= POICostUB)
			{
				// Attempt to insert the route into the skyline.
				if(LS->insertElement(ElementSkyline(el.first, POIcost), el.second.elSeq))
				{
					POICostUB = POIcost;

					// DEBUG.
//					std::cout << "Found non dominated route! ";
//					for(auto elem : el.second.elSeq)
//						std::cout << elem << " ";
//					std::cout << "(W " << el.first << ", C "<< POIcost << ")\n";
				}
			}

			// Terminate the algorithm when we find the first route
			// with the minimum possible travel cost.
			/*if(POIcost == minCost)
			{
				// ubTD = el.first;
				// std::cout << "EARLY\n";
				break;
			}*/
		}
		// 2 - Otherwise, if the sequence is below the POI cost upper bound,
		//     extend it with the NN of the last COI in the sequence.
		else
		{
			// if(POIcost <= POICostUB)
			{
				const auto NN = (el.second.elSeq.size() - 1) == 0 ?
									  this->INELookup(el.second.elSeq.back(),	// Source node.
									  el.second.elSeq.size() - 1,   	// Target category.
									  el.second.elK.back()) :			// Find k-th NN.
									  this->INEPrecomp(el.second.elSeq.back(),		// Source node.
									  el.second.elSeq.size() - 1,   			    // Target category.
									  el.second.elK.back());					    // Find k-th NN.


				// 3 - Generate the new partial sequence;
				Sequence tmpSeq = el.second;
				tmpSeq.elSeq.push_back(NN.first); tmpSeq.elW.push_back(tmpSeq.elW.back() + NN.second); tmpSeq.elK.push_back(0);
				queuePaths.push(ElQueue(el.first + NN.second, tmpSeq));

				// DEBUG.
//							std::cout << "Computing NN\n";
				//			std::cout << "Generated sequence: ";
				//			for(auto elem : tmpSeq.elSeq)
				//				std::cout << elem << " ";
				//			std::cout << "(W " << el.first + NN.second << ")\n";
			}
//			else
//				std::cout << "Not extending popped sequence\n";
		}


		// 4 - Look for the next NN for the predecessor of el. Then, add it to the queue.
		if(el.second.elSeq.size() > 1)
		{
			Sequence& seq = el.second;

			seq.elSeq.pop_back(); seq.elK.pop_back(); seq.elW.pop_back(); // Remove the last category from the sequence.
			seq.elK.back()++; 		   // Update the k of the next NN to find.
			el.first = seq.elW.back(); // Update the partial weight in el accordingly.

			double POIcostPred = 0;
			for(uint32_t i = 1; i < seq.elSeq.size(); i++)
				POIcostPred += this->getCostPOI(i-1, seq.elSeq[i]);
//			std::cout << "Attempt to generate prec. seq.: " << POIcostPred << "\n";

			// Find the k-th NN (if possible).
//			const auto kappa = seq.elK.back();
//			std::cout << "TEST: " << seq.elK.size() - 1 << "," << this->sortedVecCostPOIs[seq.elK.size() - 1].size() << "\n";
			if(/*POIcostPred <= POICostUB &&*/ seq.elK.back() < this->sortedVecCostPOIs[seq.elK.size() - 1].size())
			{
//				std::cout << "Computing NN2\n";
				const auto NNPred = (seq.elSeq.size() - 1 == 0) ?
										this->INELookup(seq.elSeq.back(),	// Source node.
										seq.elSeq.size() - 1,  	    	    // Target category.
										seq.elK.back()) :					// Find k-th NN.
										this->INEPrecomp(seq.elSeq.back(),	// Source node.
										seq.elSeq.size() - 1,  	    		// Target category.
										seq.elK.back());					// Find k-th NN.


				// Update el and push it in the queue.
				seq.elSeq.push_back(NNPred.first); seq.elW.push_back(seq.elW.back() + NNPred.second); seq.elK.push_back(0);
				el.first = seq.elW.back();
				queuePaths.push(el);

				// DEBUG.
//				std::cout << "Generated prec. sequence (k=" << kappa << "): ";
//				for(auto elem : el.second.elSeq)
//					std::cout << elem << " ";
//				std::cout << "(W " << el.first << ")\n";
			}
//			else
//				std::cout << "Reached END!!\n";
		}
	}


	std::cout << "Number of elements in the skyline: " << LS->getLS().size() << "\n";
	for(auto els : LS->getLS())
		std::cout << els.first << "," << els.second << "\n";
	for(auto els : LS->getLSS())
	{
		for(auto elp : els)
			std::cout << elp << " ";
		std::cout << "\n";
	}
}

/**
 * @brief Internal test method.
 */
void UndirectedRN::testPNE()
{
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_int_distribution<> dis(0, this->vecNodes.size() - 1);


	for(uint32_t cntTest = 0; cntTest < 1000; cntTest++)
	{
		const auto idSource = dis(gen);
		std::cout << "OSR from ID " << idSource << ": ";
		auto seq = this->PNE(idSource); seq.first.erase(seq.first.begin());
		for(auto el : seq.first)
			std::cout << el << " ";
		std::cout << "(W " << seq.second << ", C " << this->getCostSequence(seq.first) << ")\n";


		double cost = this->SPDijkstra(idSource, seq.first[0]);
		for(uint32_t i = 0; i < seq.first.size(); i++)
		{
			double costPath = this->SPDijkstra(seq.first[i], seq.first[i+1]);
			cost += costPath;
	//		std::cout << "Weight " << seq.first[i] << "," << seq.first[i+1] << ": " << costPath << "\n";
		}

		if(cost != seq.second)
		{
			std::cout << "Error PNE (check Dijkstra)!\n";
			exit(1);
		}
	}
}
