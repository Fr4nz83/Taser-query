/**
 * @brief Entry point application.
 * @author Francesco Lettich
 */


// *** INCLUDES *** //

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <queue>
#include <vector>
#include <chrono>
#include <iomanip>
#include <random>

#include "boost/program_options.hpp"
#include "GeographicLib/Geodesic.hpp"
#include "lemon/smart_graph.h"

#include "LinearSkyline.hpp"
#include "SimpleSkyline.hpp"
#include "SimpleConventionalSkyline.hpp"
#include "UndirectedRN.hpp"
#include "URN.hpp"
#include "EstimateCalculator.hpp"
#include "CASR.hpp"


// *** GLOBAL CHRONO VARIABLES *** //

double total_query_time = 0;
double OSR_time = 0;
double approach_time = 0;



// *** MAIN *** //

/**
 * @brief Entry point application.
 */
int main(int argc, char* argv[])
{
	// Local main variables.
	int startVertex = -1;
	int selectApproach = 1;
	uint32_t numRounds = 1;
	uint32_t lenSequence = 0;
	uint32_t resetHashTables = 1;
	uint32_t numPOIsPerCOI = 0;
	std::string nameFileVertices;
	std::string nameFileEdges;
	std::vector<std::string> vecNameFilesPOIs;
	std::string nameFilePrecomputeTD;

	// Command line management.
	namespace po = boost::program_options;
	po::options_description desc("Options accepted by the application");
	desc.add_options()
		("vertices,v", po::value<std::string>(&nameFileVertices)->required(), "File containing the vertices of the road network")
		("edges,e", po::value<std::string>(&nameFileEdges)->required(), "File containing the edges of the road network")
		("poi,p", po::value<std::vector<std::string>>(&vecNameFilesPOIs)->required(), "A file containing the POIs of a specific COI")
		("approach,a", po::value<int>(&selectApproach)->default_value(1), "Approach used to compute a CSR query")
		("rounds,r", po::value<uint32_t>(&numRounds)->default_value(1), "Number of queries to be computed")
		("mlen,m", po::value<uint32_t>(&lenSequence)->default_value(0), "Length of the sequence associated with CSR queries")
		("ppc", po::value<uint32_t>(&numPOIsPerCOI)->default_value(0), "Maximum number of POIs per COI")
		("hash", po::value<uint32_t>(&resetHashTables)->default_value(1), "Flag that resets the hash tables each time a query is computed (0 or 1)")
		("precomputeTD", po::value<std::string>(&nameFilePrecomputeTD), "Precomputes the network distances between POIs in different COIs. Requires to specify an output file.")
		("loadTD,l", po::value<std::string>(&nameFilePrecomputeTD), "Load the file containing the precomputed travel distances between POIs of COIs. Requires to specify an input file.")
		("source,s", po::value<int>(&startVertex)->default_value(-1), "If specified, forces this vertex ID to be starting point of CSR queries.")
		("help,h", "Print help messages");

	// "Parse" the command line.
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);

	// Manage specific options.
	if(vm.count("help"))
	{
		std::cout << desc << std::endl;
		return(0);
	}

	// Raise any error if the command line contains errors.
	po::notify(vm);



	// Basic test linear skyline.
//	LinearSkyline::testLinSky();
//	SimpleConventionalSkyline::testLinSky();
//	return(1);



	// *** Instantiate main objects *** //
	// Initialize the road network and read the POIs passed through the command line.
	std::vector<const char*> filesPOIs;
	for(uint32_t i = 0; i < vecNameFilesPOIs.size(); i++) filesPOIs.push_back(vecNameFilesPOIs[i].c_str());
	// UndirectedRN rn(nameFileVertices.c_str(), nameFileEdges.c_str(), filesPOIs);
	URN rn(nameFileVertices.c_str(), nameFileEdges.c_str(), filesPOIs);
	UndirectedRN& roadNetwork = rn;
	CASR solver(roadNetwork);



	// Manage the case where we want to precompute the travel distances between POIs.
	if(vm.count("precomputeTD"))
	{
		std::cout << "PRECOMPUTE => Precomputing network distances between POIs..." << std::endl;
		std::cout << "PRECOMPUTE => Saving information in file: " << nameFilePrecomputeTD << std::endl;

		UndirectedRN rn(nameFileVertices.c_str(), nameFileEdges.c_str(), filesPOIs);
		// rn.precomputeTDPOIs(nameFilePrecomputeTD);
		rn.precomputeTDPOIsImproved(nameFilePrecomputeTD);

		return(1);
	}

	// Manage the case where we want to load the precomputed travel distances between POIs.
	if(vm.count("loadTD"))
	{
		resetHashTables = 0;
		std::cout << "INIT => Initializing hash tables with TDs between POIs..." << std::endl;
		std::cout << "PRECOMPUTE => Loading TDs from file: " << nameFilePrecomputeTD << std::endl;
		solver.loadPrecomputedTDPOIs(nameFilePrecomputeTD);
		roadNetwork.loadPrecomputedTDPOIs(nameFilePrecomputeTD);
	}



	// *** Compute the average execution time, with starting point chosen randomly *** //
	std::random_device rd; std::mt19937 gen(rd());
	std::uniform_int_distribution<> dis(0, roadNetwork.getNumVertices() - 1);
	std::cout << "INIT => Selected approach: " << selectApproach << "\n";
	std::cout << "INIT => Rounds: " << numRounds << "\n";
	std::cout << "INIT => Reset hash tables: " << (resetHashTables ? "Yes" : "No") << "\n";
	std::cout << "INIT => length sequence: " << ((lenSequence == 0) ? vecNameFilesPOIs.size() : lenSequence) << "\n";
	for(uint32_t i = 0; i < numRounds; i++)
	{
		std::cout << "\n\nCOMPUTE => Round " << i+1 << "\n";


		// If we don't want to use caching or to precompute the travel distances, reset the
		// hash tables used to take note of the travel (and Euclidean) distance between pairs of POIs.
		if(resetHashTables) solver.resetHashTables();


		// Reset the graph to the original state. Then, extract 50 POIs per each COI and set the various
		// data structures accordingly (works only with the class URN!).
		auto t1 = std::chrono::high_resolution_clock::now();
		roadNetwork.regenerateGraph(numPOIsPerCOI);
		auto t2 = std::chrono::high_resolution_clock::now();
		std::cout << "INIT => Resetting graph: " << std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count() << " sec.\n";


		uint32_t idSource = (startVertex < 0) ? dis(gen) : startVertex;
		auto start_scoring = std::chrono::high_resolution_clock::now();
		if(selectApproach > 0) // single/multi-threaded with local linear skylines.
		{
			solver.computeQueryParallel(idSource, selectApproach, lenSequence);
		}
		else if(selectApproach == 0) // Baseline based on PNE+INE.
		{
			solver.PNE(idSource, lenSequence);
		}
		else if(selectApproach == -1) // Baseline based on PNE+INE.
		{
			solver.computeBaselinePNE(idSource, lenSequence);
		}
		else if(selectApproach == -2) // Baseline based on PNE+INE + precomputed POIs.
		{
			solver.computeBaselinePNEPrecomp(idSource, lenSequence);
		}
		else if(selectApproach == -3) // Baseline based on PNE+INE + precomputed POIs that uses conventional skylines.
		{
			solver.computeBaselinePNEPrecompConvSkyline(idSource, lenSequence);
		}
		else if(selectApproach == -4) // Baseline based on PNE+INE + precomputed POIs that does not use the pruning framework.
		{
			solver.computeBaselinePNEPrecompWithoutPruning(idSource, lenSequence);
		}
		auto end_scoring = std::chrono::high_resolution_clock::now();
		total_query_time += std::chrono::duration_cast<std::chrono::duration<double>>(end_scoring - start_scoring).count();
		std::cout << "COMPUTE => Avg. app. execution time so far: " << approach_time/(i+1) << " sec.\n";
		std::cout << "COMPUTE => Avg. tot. execution time so far: " << total_query_time/(i+1) << " sec.\n";
	}
	std::cout << "Avg. query time: " << total_query_time/numRounds << " sec.\n";
	std::cout << "Avg. OSR time: " << OSR_time/numRounds << " sec.\n";
	std::cout << "Avg. approach time: " << approach_time/numRounds << " sec.\n";
	std::cerr << OSR_time/numRounds << " " << approach_time/numRounds << " " << total_query_time/numRounds;


	return(1);
}
