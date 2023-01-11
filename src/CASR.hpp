#pragma once


// *** INCLUDES *** //

#include "UndirectedRN.hpp"
#include "LinearSkyline.hpp"
#include "SimpleSkyline.hpp"
#include "tbb/concurrent_unordered_map.h"

#include <queue>
#include <mutex>
#include <atomic>
#include <chrono>


// *** GLOBAL EXTERNAL VARIABLES *** //

extern double OSR_time;
extern double approach_time;



/**
 * @brief This class implements the main algorithm of the paper.
 */
class CASR
{
protected:

	// *** PROTECTED CLASSES *** //

	// Generic comparator used for priority queues.
	template<typename K>
	struct CmpPQueue {inline bool operator()(const K& lhs, const K& rhs) const {return lhs.first >= rhs.first;}};



	// *** PROTECTED TYPEDEFS *** //

	typedef UndirectedRN::POI POI;
	typedef std::pair<double, std::vector<uint32_t>> ElQueue;
	typedef std::priority_queue<ElQueue, std::vector<ElQueue>, CmpPQueue<ElQueue>> PriorityQueueCand;



	// *** PROTECTED FIELDS *** //

	UndirectedRN& RN;

	std::pair<std::vector<uint32_t>, double> ubPC; // Sequence representing the upper bound of the POI cost.
	std::pair<std::vector<uint32_t>, double> ubTD; // Sequence representing the upper bound of the travel distance.

	std::vector<std::vector<POI>> sortedVecCostPOIs;
	tbb::concurrent_unordered_map<uint32_t, tbb::concurrent_unordered_map<uint32_t, double>> tableTDs;
	tbb::concurrent_unordered_map<uint32_t, tbb::concurrent_unordered_map<uint32_t, double>> tableEDs;

	// Data structures used by parallel threads.
	std::vector<std::mutex> mu; // Mutexes used to manage shared skylines among multiple threads.
	std::atomic<uint32_t> idxPar;



	// *** PROTECTED METHODS *** //

	// Utility methods dealing with sequences of categories.
	inline uint64_t idCombination(const std::vector<std::vector<POI>>& vecPOIs,
 	   	   	   	   				  const std::vector<uint32_t>& it);

	inline double POICostCombination(const std::vector<std::vector<POI>>& vecPOIs,
			   	   	   	   	   	   	 const std::vector<uint32_t>& it);
	inline double POICostCombination(const std::vector<std::vector<POI>>& vecPOIs,
				   	   	   	   	   	 const std::vector<uint32_t>& it,
									 const uint32_t& startIdx,
									 const uint32_t& endIdx);
	inline double EDCombination(const std::vector<std::vector<POI>>& vecPOIs,
						 	 	const uint32_t& idSource,
								const std::vector<uint32_t>& it);
	inline double EDCombination(const std::vector<std::vector<POI>>& vecPOIs,
						   	   	   const uint32_t& idSource,
								   const std::vector<uint32_t>& it,
								   const uint32_t& startIdx,
								   const uint32_t& endIdx);
	inline double TDCombination(const std::vector<std::vector<POI>>& vecPOIs,
						 	    const uint32_t& idSource,
								const std::vector<uint32_t>& it);
	inline double TDCombination(const std::vector<std::vector<POI>>& vecPOIs,
								const uint32_t& idSource,
								const std::vector<uint32_t>& it,
								const uint32_t& startIdx,
								const uint32_t& endIdx);

	// Main methods used during the computation of CASR queries.
	std::pair<std::vector<uint32_t>, double> findPOICostUpperBound(const uint32_t& idSource);
	std::pair<std::vector<uint32_t>, double> findTravelDistUpperBound(const uint32_t& idSource);
	void filterPOIs(const double& ubTD, const double& ubPC, const uint32_t& idSource);

	// Methods that implement the parallel approach.
	std::vector<uint32_t> getNextSlice(const std::vector<uint32_t>& currSlice,
									   const uint32_t& idSource,
									   LinearSkyline& skyline,
									   const uint32_t& numThreads,
									   const uint32_t& tid);
	void exploreCandidateSpaceParallel(const std::pair<std::vector<uint32_t>, double>& ubPC,
									   const std::pair<std::vector<uint32_t>, double>& ubTD,
									   const uint32_t& idSource,
									   const uint32_t& numThreads,
									   const uint32_t& tid,
									   std::vector<LinearSkyline>& LS);



public:

	// *** PUBLIC CTORS/DTOR *** //

	CASR(UndirectedRN& RN) : RN(RN) {};
	~CASR(){};



	// *** PUBLIC METHODS *** //

	// Main methods.
	void resetHashTables();
	void loadPrecomputedTDPOIs(const std::string& nameFile);

	void computeQueryParallel(const uint32_t& idSource, const uint32_t& numThreads = 1, const uint32_t& lenSequence = 0);
	void PNE(const uint32_t& idSource, const uint32_t& lenSequence = 0);
	void computeBaselinePNE(const uint32_t& idSource, const uint32_t& lenSequence = 0);
	void computeBaselinePNEPrecomp(const uint32_t& idSource, const uint32_t& lenSequence = 0);
	void computeBaselinePNEPrecompConvSkyline(const uint32_t& idSource, const uint32_t& lenSequence = 0);
	void computeBaselinePNEPrecompWithoutPruning(const uint32_t& idSource, const uint32_t& lenSequence = 0);
};



// *** INLINE PROTECTED METHODS DEFINITIONS *** //

/**
 * @brief Generate a unique ID for each candidate.
 *
 * @note Formula: id_1 + (id_2-1)*(d_1) + (id_3-1)*(d_2*d1) + (id_4-1)*(d_3*d_2*d_1),
 *  	 where id_X is the position in dimension X, while d_Y is the number of elements of
 *  	 dimension Y.
 */
inline uint64_t CASR::idCombination(const std::vector<std::vector<POI>>& vecPOIs,
		   	   	   	   	   	   	   	const std::vector<uint32_t>& it)
{
	const uint32_t numTypes = vecPOIs.size();
	uint64_t runningSum = it[0];
	uint64_t size = 1;
	for(uint32_t i = 1; i < numTypes; i++)
	{
		size *= vecPOIs[i-1].size();
		runningSum += size * it[i];
	}

	// DEBUG.
//	std::cout << "ID comb ";
//	for(uint32_t i = 0; i < numTypes; i++)
//		std::cout << it[i] << " ";
//	std::cout << ": " << runningSum + it[numTypes - 1] << "\n";

	return(runningSum);
}

inline double CASR::POICostCombination(const std::vector<std::vector<POI>>& vecPOIs,
									   const std::vector<uint32_t>& it)
{
	const double numTypes = vecPOIs.size();
	double cost = 0;
	for(uint32_t i = 0; i < numTypes; i++)
		cost += vecPOIs[i][it[i]].first;
	return(cost);
}

inline double CASR::POICostCombination(const std::vector<std::vector<POI>>& vecPOIs,
									   const std::vector<uint32_t>& it,
									   const uint32_t& startIdx,
									   const uint32_t& endIdx)
{
	double cost = 0;
	for(uint32_t i = startIdx; i < endIdx; i++)
		cost += vecPOIs[i][it[i]].first;
	return(cost);
}

/**
 * @brief This method computes the Euclidean distance of an instance of a sequence of POI categories.
 */
inline double CASR::EDCombination(const std::vector<std::vector<POI>>& vecPOIs,
						   	      const uint32_t& idSource,
								  const std::vector<uint32_t>& it)
{
	double cost = 0;
	for(uint32_t type = 0; type < vecPOIs.size(); type++)
	{
		// Get the IDs of the vertices.
		const uint32_t source = (type == 0) ? idSource : vecPOIs[type - 1][it[type - 1]].second;
		const uint32_t target = vecPOIs[type][it[type]].second;

		// Check if the TD was already computed.
		const auto lookupTarget = this->tableEDs[source].find(target);
		if(lookupTarget != this->tableEDs[source].end())
		{
//			std::cout << "Lookup! " << source << "," << target << ": " << lookupTarget->second << "\n";
			cost += lookupTarget->second;
		}
		else
		{
			const double partCost = this->RN.computeDistancePointsIDs(source, target);
			cost += partCost;

			// Update the hash-table.
			this->tableEDs[source][target] = partCost;
		}
	}

	return(cost);
}

inline double CASR::EDCombination(const std::vector<std::vector<POI>>& vecPOIs,
									 const uint32_t& idSource,
									 const std::vector<uint32_t>& it,
									 const uint32_t& startIdx,
									 const uint32_t& endIdx)
{
	double cost = 0;
	for(uint32_t type = startIdx; type < endIdx; type++)
	{
		// Get the IDs of the vertices.
		const uint32_t source = (type == 0) ? idSource : vecPOIs[type - 1][it[type - 1]].second;
		const uint32_t target = vecPOIs[type][it[type]].second;

		// Check if the ED was already computed.
		const auto lookupTarget = this->tableEDs[source].find(target);
		if(lookupTarget != this->tableEDs[source].end())
		{
//			std::cout << "Lookup! " << source << "," << target << ": " << lookupTarget->second << "\n";
			cost += lookupTarget->second;
		}
		else
		{
			const double partCost = this->RN.computeDistancePointsIDs(source, target);
			cost += partCost;

			// Update the hash-table.
			this->tableEDs[source][target] = partCost;
		}
	}

	return(cost);
}

/**
 * @brief This method computes the travel distance of an instance of a sequence of POI categories.
 */
inline double CASR::TDCombination(const std::vector<std::vector<POI>>& vecPOIs,
						   const uint32_t& idSource,
						   const std::vector<uint32_t>& it)
{
	double cost = 0;
	for(uint32_t type = 0; type < vecPOIs.size(); type++)
	{
		// Get the IDs of the vertices.
		const uint32_t source = (type == 0) ? idSource : vecPOIs[type - 1][it[type - 1]].second;
		const uint32_t target = vecPOIs[type][it[type]].second;

		// Check if the TD was already computed.
		const auto lookupTarget = this->tableTDs[source].find(target);
		if(lookupTarget != this->tableTDs[source].end())
		{
//			std::cout << "Lookup! " << source << "," << target << ": " << lookupTarget->second << "\n";
			cost += lookupTarget->second;
		}
		else
		{
			const double partCost = this->RN.SPDijkstra(source, target);
			cost += partCost;

			// Update the hash-table.
			this->tableTDs[source][target] = partCost;

			// std::cout << "Dijkstra1! " << source << "," << target << ": " << partCost << "\n";
		}
	}

	return(cost);
}

/**
 * @brief This method computes the travel distance of an instance of a sequence of POI categories.
 */
inline double CASR::TDCombination(const std::vector<std::vector<POI>>& vecPOIs,
								  const uint32_t& idSource,
								  const std::vector<uint32_t>& it,
								  const uint32_t& startIdx,
								  const uint32_t& endIdx)
{
	double cost = 0;
	for(uint32_t type = startIdx; type < endIdx; type++)
	{
		// Get the IDs of the vertices.
		const uint32_t source = (type == 0) ? idSource : vecPOIs[type - 1][it[type - 1]].second;
		const uint32_t target = vecPOIs[type][it[type]].second;

		// Check if the TD was already computed.
		const auto lookupTarget = this->tableTDs[source].find(target);
		if(lookupTarget != this->tableTDs[source].end())
		{
//			std::cout << "Lookup! " << source << "," << target << ": " << lookupTarget->second << "\n";
			cost += lookupTarget->second;
		}
		else
		{
			const double partCost = this->RN.SPDijkstra(source, target);
			cost += partCost;

			// Update the hash-table.
			this->tableTDs[source][target] = partCost;

			// std::cout << "Dijkstra1! " << source << "," << target << ": " << partCost << "\n";
		}
	}

	return(cost);
}
