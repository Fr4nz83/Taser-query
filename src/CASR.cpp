// *** INCLUDES *** //

#include "CASR.hpp"

#include <list>
#include <thread>
#include <iostream>
#include <fstream>


// *** PROTECTED METHODS DEFINITIONS *** //

/**
 * @brief This method eliminates those POIs that are surely dominated by the upper bounds.
 */
void CASR::filterPOIs(const double& ubTD, const double& ubPC, const uint32_t& idSource)
{
	std::cout << "Filtering with ubTD " << ubTD << " and ubPC " << ubPC << " (ID " << idSource << ")\n";
	for(uint32_t i = 0; i < this->sortedVecCostPOIs.size(); i++)
	{
		std::cout << "Filtering category " << i << " (" << this->sortedVecCostPOIs[i].size() << " POIs)\n";
		// Accumulate the POIs in an ordered set.
		std::set<POI> tmpSet;
		for(auto el : this->sortedVecCostPOIs[i])
		{
			const double Edist = this->RN.computeDistancePointsIDs(idSource, el.second);
//			if(el.first <= ubPC && Edist <= ubTD) tmpSet.insert(el);
			if(Edist <= ubTD) tmpSet.insert(el);
//			else std::cout << el.second << " rejected!\n";
		}

		// Update the vector of sorted POIs for this category.
		this->sortedVecCostPOIs[i] = std::vector<POI>(tmpSet.begin(), tmpSet.end());
		std::cout << "# non dominated POIs category " << i << ": " << this->sortedVecCostPOIs[i].size() << "\n";
	}
}

/**
 * @brief This method solves the OSR query, i.e., finds the sequence yielding the minimum
 * 		  travel distance (and thus the POI cost upper bound).
 *
 * @return The OSR.
 */
std::pair<std::vector<uint32_t>, double> CASR::findPOICostUpperBound(const uint32_t& idSource)
{
	// Compute the OSR.
	auto seq = this->RN.PNE(idSource);

	// Delete the ID of the source from the sequence.
	seq.first.erase(seq.first.begin());

	// Return the OSR.
	return(seq);
}

/**
 * @brief This method computes the path yielding the minimum POI cost, and thus the upper bound
 * 		  for the travel distance.
 */
std::pair<std::vector<uint32_t>, double> CASR::findTravelDistUpperBound(const uint32_t& idSource)
{
	// 1 - Get the candidate route having minimum absolute cost and compute its travel distance.
	std::vector<uint32_t> itCurr(this->RN.getNumTypesPOIs(), 0);
//	double minDist = this->TDCombination(this->sortedVecCostPOIs,
//										  idSource,
//										  itCurr);

	// Replace the indexes within the candidate with the IDs of the vertices they refer to.
	for(uint32_t type = 0; type < itCurr.size(); type++)
		itCurr[type] = this->sortedVecCostPOIs[type][itCurr[type]].second;

	// Return the candidate and its travel distance.
	return(std::pair<std::vector<uint32_t>, double>(itCurr, 0));
}

/**
 * @brief This method explore the candidate space, looking for linearly non-dominated paths.
 *
 * @param ubPC upper bound POI cost.
 * @param ubTD upper bound travel distance.
 * @param source Query's source.
 */
void CASR::exploreCandidateSpaceParallel(const std::pair<std::vector<uint32_t>, double>& ubPC,
									     const std::pair<std::vector<uint32_t>, double>& ubTD,
									     const uint32_t& idSource,
										 const uint32_t& numThreads,
										 const uint32_t& tid,
										 std::vector<LinearSkyline>& LS)
{
	// TODO: this method doesn't work anymore as "ubTravelDist" is not anymore computed before
	//		 executing this part!
	// 1 - Initialize the data structures needed to search the candidate space.
	const double ubPOICost = this->RN.getCostSequence(ubPC.first);
	const double ubTravelDist = this->ubTD.second;
	const uint32_t numTypes = this->sortedVecCostPOIs.size();
	const uint32_t idxLastType = numTypes - 1;


	// *** Skyline initialization *** //
	// Manage the particular case where the travel cost of the two upper bounds is the same:
	// we can immediately return a skyline containing only the route determining the travel cost upper
	// as it is the one with minimum travel distance (and hence the solution).
	const LinearSkyline::ElementSkyline upperTD(this->RN.getCostSequence(ubTD.first), ubTD.second);
	const LinearSkyline::ElementSkyline upperPC(this->RN.getCostSequence(ubPC.first), ubPC.second);
	if(upperTD.first == upperPC.first)
	{
		LS[tid] = LinearSkyline(upperPC, ubPC.first, upperPC, ubPC.first);
		return;
	}
	this->mu[tid].lock();
	LS[tid] = LinearSkyline(upperTD, ubTD.first, upperPC, ubPC.first);
	LinearSkyline LSLocal = LS[tid];
	this->mu[tid].unlock();


	// These variables regulate the frequency with which threads share skylines.
	// const uint32_t blockCandidates = 200;
	// const uint32_t blockSlices = std::max((uint32_t) 1, blockCandidates / (uint32_t) this->sortedVecCostPOIs[idxLastType].size());
	const uint32_t blockSlices = 64;


	if(tid == 0) std::cout << "Debug: block slices: " << blockSlices << "\n";
	if(tid == 0) std::cout << "Debug: travel dist ub: " << ubTravelDist << "\n";
	if(tid == 0) std::cout << "Debug: POI cost ub: " << ubPOICost << "\n";
//	std::cout << "Debug: TID " << tid << "\n";
//	uint32_t numCandidates = 0;
	uint32_t cntSlices = 0;
	std::vector<uint32_t> vecIndexes(numTypes, 0); vecIndexes[0] = this->idxPar++;
	// Prevent the execution of threads that go outside the range of families of slices available.
	if(vecIndexes[0] >= this->sortedVecCostPOIs[0].size()) vecIndexes.clear();
	while(vecIndexes.size())
	{
		// DEBUG -- verify that slices are generated in lexicographic order.
//		std::cout << "Current sequence: ";
//		for(auto index : vecIndexes)
//			std::cout << index << " ";
//		std::cout << "\n";


		// Verify that this slice can generate non dominated routes.
		// To this end, check first the Euclidean distance, and then the travel distance.
		// Check also the POI cost.
		double POICostCombination = this->POICostCombination(this->sortedVecCostPOIs, vecIndexes);
		const double EDCombination = this->EDCombination(this->sortedVecCostPOIs, idSource, vecIndexes, 0, idxLastType);

		// Determine the initial travel distance upper bound for this slice.
		auto ptr = LSLocal.findRightNeighbor(LinearSkyline::ElementSkyline(POICostCombination, 0));
		double ubTDLocal = (ptr == LSLocal.cbegin()) ? ptr->second : (--ptr)->second;

		// DEBUG.
//			std::cout << "Evaluating slice! CP:" << POICostCombination << ", ED:" << EDCombination << "\n";
//			std::cout << "ubTDLocal: " << ubTDLocal << "\n";

		POICostCombination -= this->sortedVecCostPOIs[idxLastType][vecIndexes[idxLastType]].first;
		while(vecIndexes[idxLastType] < this->sortedVecCostPOIs[idxLastType].size())
		{
			const double POICost = POICostCombination + this->sortedVecCostPOIs[idxLastType][vecIndexes[idxLastType]].first;
			const double Edist = EDCombination + this->EDCombination(this->sortedVecCostPOIs,
																	 idSource,
																	 vecIndexes,
																	 idxLastType,
																	 idxLastType+1);

			// Check if we can early terminate the evaluation of the current slice.
			if(POICost > ubPOICost) break;

			// Evaluate the current candidate.
			if(Edist <= ubTDLocal)
			{
				// Compute the actual travel distance of this candidate.
				const double dist = this->TDCombination(this->sortedVecCostPOIs,
														idSource,
														vecIndexes);

				// Since we are checking if this element can be inserted into the skyline,
				// we can spend a few operations to update the travel distance upper bound
				// from the skyline (gives some performance advantages).
				ptr = LSLocal.findRightNeighbor(LinearSkyline::ElementSkyline(POICost, 0));
				ubTDLocal = (ptr == LSLocal.cbegin()) ? ptr->second : (--ptr)->second;

				// Attempt to insert the current candidate in the local skyline.
				// If the element is inserted, update the travel distance upper bound.
				if(LSLocal.insertElement(LinearSkyline::ElementSkyline(POICost, dist), vecIndexes))
					ubTDLocal = dist;
			}

			// Get the next candidate in this slice.
			vecIndexes.back()++;

			// DEBUG -- print the number of candidates evaluated so far.
//				numCandidates++;
//				if(numCandidates % 100 == 0) std::cout << "Candidates eval: " << numCandidates << "\n";
		}


		// *** BEGIN => Periodically update the "shared" skylines between the threads. *** //
		if(cntSlices++ % blockSlices == 0)
		{
			for(uint32_t i = 0; i < numThreads; i++)
			{
				if(i != tid)
				{
					this->mu[i].lock();
					LinearSkyline tmp = LS[i];
					this->mu[i].unlock();

					LSLocal.mergeSkyline(tmp);
				}
			}
			this->mu[tid].lock();
			LS[tid] = LSLocal;
			this->mu[tid].unlock();
		}
		// *** END => Periodically update the "shared" skylines between the threads. *** //

		// Get the next slice to process.
		vecIndexes = getNextSlice(vecIndexes, idSource, LSLocal, numThreads, tid);
	}


	// Write the final skyline into the map.
	this->mu[tid].lock();
	LS[tid] = LSLocal;
	this->mu[tid].unlock();
	// std::cout << "Termina TID " << tid << "\n";
}

/**
 * @brief This method returns the next useful slice of candidates (if any).
 * @note The method returns an empty vector if no more slices are available for processing.
 */
std::vector<uint32_t> CASR::getNextSlice(const std::vector<uint32_t>& currSlice,
										 const uint32_t& idSource,
										 LinearSkyline& skyline,
										 const uint32_t& numThreads,
										 const uint32_t& tid)
{
	// Particular case: if the sequence of COIs has length 1, we do not need to find any other slice.
	// Return an empty vector.
	if(currSlice.size() == 1) return(std::vector<uint32_t>(0));

	// Determine the new slice: reset the last index and start from the penultimate
	// element in the sequence.
	std::vector<uint32_t> newSlice(currSlice); newSlice.back() = 0;
	const uint32_t penultimate = newSlice.size() - 2;
	uint32_t indexElSeq = penultimate;
	bool done = false;
	while(!done)
	{
		// Update the counter of the current index.
		newSlice[indexElSeq] = (indexElSeq != 0) ? newSlice[indexElSeq] + 1 : this->idxPar++;

		// Case 1: we are within the range of the current index.
		if(newSlice[indexElSeq] < this->sortedVecCostPOIs[indexElSeq].size())
		{
			// Flag used to terminate the cycle.
			done = true;

			// Calculate the POI Cost and the Euclidean distance associated with
			// the subsequence between idx "0" and idx "indexElSeq".
			const double POICostCombination = this->POICostCombination(this->sortedVecCostPOIs, newSlice);
			const double EDCombination = this->EDCombination(this->sortedVecCostPOIs,
															 idSource, newSlice,
															 0, indexElSeq+1);

			const double ubPC = skyline.getSecondUpperBound().first;


			// Case 1: POI cost greater than upper bound.
			if(POICostCombination > ubPC)
			{
				// Case where it is not possible to find more slices: return an empty vector.
				if(indexElSeq == 0) return(std::vector<uint32_t>(0));

				newSlice[indexElSeq] = 0;
				indexElSeq--;
				done = false;
			}
			// Case 2: Euclidean distance greater than the travel upper bound.
			else
			{
				auto ptr = skyline.findRightNeighbor(LinearSkyline::ElementSkyline(POICostCombination, 0));
				const double ubTD = (ptr == skyline.cbegin()) ? ptr->second : (--ptr)->second;
				if(EDCombination > ubTD)
				{
					done = false;
				}
			}
		}
		// Case 2: we are outside the range of the current index. Reset it to 0 and go
		// to the preceding index (if possible).
		else
		{
			// Case where it is not possible to find more slices: return an empty vector.
			if(indexElSeq == 0) return(std::vector<uint32_t>(0));

			newSlice[indexElSeq] = 0;
			indexElSeq--;
		}
	}

	// Return the new slice.
	return(newSlice);
}



// *** PUBLIC METHODS DEFINITIONS *** //

/**
 * @brief Method used to reset the lookup tables used to speedup the computation of Euclidean distances
 * 		  and travel distances.
 */
void CASR::resetHashTables()
{
	this->tableEDs.clear();
	this->tableTDs.clear();
}

void CASR::loadPrecomputedTDPOIs(const std::string& nameFile)
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
		this->tableEDs[source][target] = dist;
		this->tableEDs[target][source] = dist;
		this->tableTDs[source][target] = dist;
		this->tableTDs[target][source] = dist;
	}

	infile.close();
}

/**
 * @brief Method that computes a CSR query via a parallel approach that does not rely on a priority queue.
 */
void CASR::computeQueryParallel(const uint32_t& idSource, const uint32_t& numThreads, const uint32_t& lenSequence)
{
	std::cout << "COMPUTE => ID starting vertex " << idSource << "\n";


	// 1 - Copy the vectors of sorted POIs from RN.
	std::cout << "COMPUTE => Initialize data structures...\n";
	this->RN.setupQuerySequence(lenSequence);
	this->sortedVecCostPOIs = this->RN.getVecSortedPOIs();


	// 2 - Compute the POI cost upper bound.
	std::cout << "COMPUTE => Finding upper bound POI cost (minimum travel distance)...\n";

	auto start_scoring = std::chrono::high_resolution_clock::now();
	this->ubPC = this->findPOICostUpperBound(idSource);
	auto end_scoring = std::chrono::high_resolution_clock::now();
	OSR_time += std::chrono::duration_cast<std::chrono::duration<double>>(end_scoring - start_scoring).count();

	for(auto el : this->ubPC.first)
		std::cout << el << " ";
	std::cout << "(W " << this->ubPC.second << ", C " << this->RN.getCostSequence(this->ubPC.first) << ")\n";
	const double ubPOICost = this->RN.getCostSequence(this->ubPC.first);


	// 3 - Compute the travel distance upper bound.
	std::cout << "COMPUTE => Finding upper bound travel distance (minimum POI Cost)...\n";
	this->ubTD = this->findTravelDistUpperBound(idSource);
	for(auto el : this->ubTD.first)
		std::cout << el << " ";
	std::cout << "(W " << this->ubTD.second << ", C " << this->RN.getCostSequence(this->ubTD.first) << ")\n";
	const double ubTravelDist = this->ubTD.second;


	// 4 - Filter out from the local copy of the vector of sorted POIs those POIs which Euclidean distance
	//	   from idSource is above the length of the TD upper bound (they are surely dominated), and those
	//	   that have a cost greater than the PC upper bound.
	std::cout << "COMPUTE => Pruning POIs that surely yield dominated paths...\n";
	this->filterPOIs(ubTravelDist, ubPOICost, idSource);


	// 5 - Explore the candidate search space.
	std::cout << "COMPUTE => Linear skyline...\n";
	start_scoring = std::chrono::high_resolution_clock::now();


	// 6 - Allocate threads' resources.
	const uint32_t maxThreads = (this->sortedVecCostPOIs[0].size() > numThreads) ? numThreads : this->sortedVecCostPOIs[0].size();
	std::vector<LinearSkyline> LS(maxThreads);
	std::vector<std::thread> t(maxThreads);
	this->mu = std::vector<std::mutex>(maxThreads);
	this->idxPar = 0;


	// 7 - Start the threads.
	for(uint8_t tid = 0; tid < maxThreads; tid++)
		t[tid] = std::thread(&CASR::exploreCandidateSpaceParallel,
							 this,
							 std::ref(this->ubPC),
							 std::ref(this->ubTD),
							 std::ref(idSource),
							 std::ref(maxThreads),
							 tid,
							 std::ref(LS));
	// All the threads have to reach the join below, before allowing the application to continue with
	// its execution.
	for(uint8_t tid = 0; tid < maxThreads; tid++) t[tid].join();


	// DEBUG -- print the current linear skyline.
//	for(uint32_t i = 0; i < LS.size(); i ++)
//	{
//		std::cout << "RESULT => Number of elements in the " << i+1 << "-th skyline: " << LS[i].getLS().size() << "\n";
//		for(auto el : LS[i].getLS())
//			std::cout << el.first << "," << el.second << "\n";
//		std::cout << "\n";
//	}


	// Merge the skylines.
	std::cout << "RESULT => Merging " << maxThreads << " skylines...\n";
	for(uint32_t i = 1; i < LS.size(); i++) LS[0].mergeSkyline(LS[i]);
	end_scoring = std::chrono::high_resolution_clock::now();
	approach_time += std::chrono::duration_cast<std::chrono::duration<double>>(end_scoring - start_scoring).count();


	// DEBUG -- print the final linear skyline.
	std::cout << "RESULT => Number of elements in the final skyline: " << LS[0].getLS().size() << "\n";
	for(auto el : LS[0].getLS())
	{
		std::cout << el.first << "," << el.second;
		std::cout << " -- Sequence: ";
		auto seq = LS[0].getSequence(el.first);
		for(auto el : seq) std::cout << el << " ";
		std::cout << "\n";
	}
	std::cout << "\n";
}

void CASR::PNE(const uint32_t& idSource, const uint32_t& lenSequence)
{
	// 1 - Copy the vectors of sorted POIs from RN.
	std::cout << "COMPUTE => Initialize data structures...\n";
	this->RN.setupQuerySequence(lenSequence);
	this->sortedVecCostPOIs = this->RN.getVecSortedPOIs();


	// 2 - Initialize the hash table where each POI contains the k NN with respect
	//	   to each type of POI used in the query.
	std::cout << "COMPUTE => Initializing lookup table k nearest neighbors POIs...\n";
	this->RN.resetHashTableTDs();


	// 3 - Compute the kNNs of the source.
	std::cout << "COMPUTE => Computing k-nearest neighbors source...\n";
	this->RN.initkNNsSource(idSource);


	// 4 - Execute PNE.
	std::cout << "COMPUTE => Executing PNE.\n";
	auto res = this->RN.PNE(idSource);


	std::cout << "OSR: ";
	for(const auto& el : res.first)
		std::cout << el << " ";
	std::cout << "- Length: " << res.second << "\n";
}

void CASR::computeBaselinePNE(const uint32_t& idSource, const uint32_t& lenSequence)
{
	// 1 - Copy the vectors of sorted POIs from RN.
	std::cout << "COMPUTE => Initialize data structures...\n";
	this->RN.setupQuerySequence(lenSequence);
	this->sortedVecCostPOIs = this->RN.getVecSortedPOIs();


	// 2 - Initialize the hash table where each POI contains the k NN with respect
	//	   to each type of POI used in the query.
	std::cout << "COMPUTE => Initializing lookup table k nearest neighbors POIs...\n";
	this->RN.resetHashTableTDs();


	// Compute the kNNs of the source.
	std::cout << "COMPUTE => Computing k-nearest neighbors source...\n";
	this->RN.initkNNsSource(idSource);


	// 2 - Compute the travel distance upper bound.
	std::cout << "COMPUTE => Finding upper bound travel distance (minimum travel Cost)...\n";
	auto start_scoring = std::chrono::high_resolution_clock::now();
	this->ubTD = this->findTravelDistUpperBound(idSource);
	auto end_scoring = std::chrono::high_resolution_clock::now();
	OSR_time += std::chrono::duration_cast<std::chrono::duration<double>>(end_scoring - start_scoring).count();
	for(auto el : this->ubTD.first)
		std::cout << el << " ";
	std::cout << "(W " << this->ubTD.second << ", C " << this->RN.getCostSequence(this->ubTD.first) << ")\n";


	// 3 - Execute the baseline based on PNE+INE.
	this->RN.PNEBaseline(idSource, this->RN.getCostSequence(this->ubTD.first));
}

/**
 * @brief CSR_{PNE} using precomputed network distances between POIs.
 */
void CASR::computeBaselinePNEPrecomp(const uint32_t& idSource, const uint32_t& lenSequence)
{
	// 1 - Copy the vectors of sorted POIs from RN.
	std::cout << "COMPUTE => Initialize data structures...\n";
	this->RN.setupQuerySequence(lenSequence);
	this->sortedVecCostPOIs = this->RN.getVecSortedPOIs();


	// 2 - Reset the table containing the POIs' kNNs.
	std::cout << "COMPUTE => Resetting hash table k-nearest neighbors POIs...\n";
	this->RN.resetHashTableTDs();


	// 3 - Initialize the hash table where each POI contains the k NN with respect
	//	   to each type of POI used in the query.
	std::cout << "COMPUTE => Initializing k nearest neighbors POIs...\n";
	this->RN.initPrecomputedNNsPOIs();


	// Compute the kNNs of the source.
	std::cout << "COMPUTE => Computing k-nearest neighbors source...\n";
	auto start_scoring = std::chrono::high_resolution_clock::now();
	this->RN.initkNNsSource(idSource);
	auto end_scoring = std::chrono::high_resolution_clock::now();
	OSR_time += std::chrono::duration_cast<std::chrono::duration<double>>(end_scoring - start_scoring).count();


	// 4 - Compute the travel distance upper bound.
	std::cout << "COMPUTE => Finding upper bound travel distance (minimum POI Cost)...\n";
	this->ubTD = this->findTravelDistUpperBound(idSource);
	for(auto el : this->ubTD.first)
		std::cout << el << " ";
	std::cout << "(W " << this->ubTD.second << ", C " << this->RN.getCostSequence(this->ubTD.first) << ")\n";


	// 5 - Execute the approach.
	this->RN.PNEBaselinePrecomp(idSource, this->RN.getCostSequence(this->ubTD.first));
}

/**
 * @brief CSR_{PNE} using precomputed network distances between POIs.
 */
void CASR::computeBaselinePNEPrecompConvSkyline(const uint32_t& idSource, const uint32_t& lenSequence)
{
	std::cout << "COMPUTE => Using conventional skylines...\n";


	// 1 - Copy the vectors of sorted POIs from RN.
	std::cout << "COMPUTE => Initialize data structures...\n";
	this->RN.setupQuerySequence(lenSequence);
	this->sortedVecCostPOIs = this->RN.getVecSortedPOIs();


	// 2 - Reset the table containing the POIs' kNNs.
	std::cout << "COMPUTE => Resetting hash table k-nearest neighbors POIs...\n";
	this->RN.resetHashTableTDs();


	// 3 - Initialize the hash table where each POI contains the k NN with respect
	//	   to each type of POI used in the query.
	std::cout << "COMPUTE => Initializing k nearest neighbors POIs...\n";
	this->RN.initPrecomputedNNsPOIs();


	// Compute the kNNs of the source.
	std::cout << "COMPUTE => Computing k-nearest neighbors source...\n";
	auto start_scoring = std::chrono::high_resolution_clock::now();
	this->RN.initkNNsSource(idSource);
	auto end_scoring = std::chrono::high_resolution_clock::now();
	OSR_time += std::chrono::duration_cast<std::chrono::duration<double>>(end_scoring - start_scoring).count();


	// 4 - Compute the travel distance upper bound.
	std::cout << "COMPUTE => Finding upper bound travel distance (minimum POI Cost)...\n";
	this->ubTD = this->findTravelDistUpperBound(idSource);
	for(auto el : this->ubTD.first)
		std::cout << el << " ";
	std::cout << "(W " << this->ubTD.second << ", C " << this->RN.getCostSequence(this->ubTD.first) << ")\n";


	// 5 - Execute the approach.
	this->RN.PNEBaselinePrecompConvSkyline(idSource, this->RN.getCostSequence(this->ubTD.first));
}

/**
 * @brief CSR_{PNE} using precomputed network distances between POIs.
 */
void CASR::computeBaselinePNEPrecompWithoutPruning(const uint32_t& idSource, const uint32_t& lenSequence)
{
	std::cout << "COMPUTE => leaving out the pruning framework...\n";


	// 1 - Copy the vectors of sorted POIs from RN.
	std::cout << "COMPUTE => Initialize data structures...\n";
	this->RN.setupQuerySequence(lenSequence);
	this->sortedVecCostPOIs = this->RN.getVecSortedPOIs();


	// 2 - Reset the table containing the POIs' kNNs.
	std::cout << "COMPUTE => Resetting hash table k-nearest neighbors POIs...\n";
	this->RN.resetHashTableTDs();


	// 3 - Initialize the hash table where each POI contains the k NN with respect
	//	   to each type of POI used in the query.
	std::cout << "COMPUTE => Initializing k nearest neighbors POIs...\n";
	this->RN.initPrecomputedNNsPOIs();


	// Compute the kNNs of the source.
	std::cout << "COMPUTE => Computing k-nearest neighbors source...\n";
	auto start_scoring = std::chrono::high_resolution_clock::now();
	this->RN.initkNNsSource(idSource);
	auto end_scoring = std::chrono::high_resolution_clock::now();
	OSR_time += std::chrono::duration_cast<std::chrono::duration<double>>(end_scoring - start_scoring).count();


	// 4 - Compute the travel distance upper bound.
	std::cout << "COMPUTE => Finding upper bound travel distance (minimum POI Cost)...\n";
	this->ubTD = this->findTravelDistUpperBound(idSource);
	for(auto el : this->ubTD.first)
		std::cout << el << " ";
	std::cout << "(W " << this->ubTD.second << ", C " << this->RN.getCostSequence(this->ubTD.first) << ")\n";


	// 5 - Execute the approach.
	this->RN.PNEBaselinePrecompWithoutPruning(idSource, this->RN.getCostSequence(this->ubTD.first));
}
