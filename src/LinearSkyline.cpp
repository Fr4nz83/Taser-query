/**
 * @brief Class modeling a linear skyline.
 * @author Francesco Lettich
 */

// *** INCLUDES *** //

#include "LinearSkyline.hpp"

#include <vector>
#include <ctime>
#include <iostream>
#include <cstdlib>
#include <fstream>



// *** PUBLIC CTORS DEFINITIONS *** //

/**
 * @brief Initializes a new linear skyline with two lower bounds: lb1, which is associated
 *  	  with the first type of cost, and lb2 which is associated with the second type of cost.
 */
LinearSkyline::LinearSkyline(const ElementSkyline& ub1, const std::vector<uint32_t>& ub1Seq,
							 const ElementSkyline& ub2, const std::vector<uint32_t>& ub2Seq)
{
	this->LS.insert(ub1); this->tableSequences[ub1.first] = ub1Seq;
	this->LS.insert(ub2); this->tableSequences[ub2.first] = ub2Seq;
}



// *** PUBLIC METHODS DEFINITIONS *** //

/**
 * @brief This method attempts to insert "el" into the linear skyline.
 *
 * @return true if the insertion is successful, false otherwise.
 * @note The method may remove elements from the linear skyline if "el" is inserted into the skyline.
 */
bool LinearSkyline::insertElement(const ElementSkyline& el, const std::vector<uint32_t>& elSeq)
{
	Skyline::const_iterator ptr = this->findRightNeighbor(el);
	if(ptr == this->LS.end()) return(false); // Case where el is greater or equal than the upper bound
											 // in the second dimension: we don't insert anything.

	const auto rneigh = *ptr--;
	const auto lneigh = *ptr;
//	std::cout << "Left neighbor: " << lneigh.first << "," << lneigh.second << "\n";
//	std::cout << "Right neighbor: " << rneigh.first << "," << rneigh.second << "\n";

	// First, check if the left neighbor conventionally dominates "el": if so, refuse the insertion.
	// NOTE: we discard also those routes that have the same travel cost and travel distance w.r.t. their
	//	     left neighbor.
//	if(this->checkDomination(lneigh, el)) return(false);
	if(lneigh.second <= el.second) return(false);

//	std::cout << "Non convenzionalmente dominato.\n";

	// Check for linear domination of "el" by its left and right neighbors:
	// if dominated, refuse the insertion.
	if(this->checkNormalDomination(lneigh, rneigh, el)) return(false);
//	std::cout << "Non linearmente dominato.\n";


	// If arrived here, insert "el" in the skyline...
	auto test = this->LS.insert(el);

	// Case where there's already an element with the same key.
	if(!test.second)
	{
		// Case 1: the element already present in the skyline dominates (or is equivalent to) the new one.
		// if(test.first->second <= el.second) return(false);

		// Case 2: the new element dominates the one in the skyline.
		this->LS.erase(test.first);
		test = this->LS.insert(el);
	}
	this->tableSequences[el.first] = elSeq;


	// Verify if there are skyline elements on the left of "el" that are linearly dominated.
//	std::cout << "Verifica se ci sono elementi da rimuovere a SX nella skyline.\n";
	auto prevPtr = test.first; prevPtr--;
	while(prevPtr != this->LS.begin())
	{
		auto prevPrevPtr = prevPtr; prevPrevPtr--;
		if(this->checkNormalDomination(*prevPrevPtr, el, *prevPtr))
		{
			// NOTE: after deletion, prevPtr points to the element coming next to the deleted one.
			this->tableSequences.erase(prevPtr->first);
			prevPtr = this->LS.erase(prevPtr); prevPtr--;
//			std::cout << "Rimozione elemento SX skyline!\n";
		}
		else
			break;
	}


	// Verify if there are skyline elements on the right of "el" that are linearly dominated.
	//	std::cout << "Verifica se ci sono elementi da rimuovere a DX nella skyline.\n";
	// NOTE: This part is not needed in case we attempt to insert elements in increasing order of
	//		 the cost associated with the first dimension.
	auto succPtr = test.first; succPtr++;
	while(succPtr != std::prev(this->LS.end()))
	{
		auto succSuccPtr = succPtr; succSuccPtr++;
		if(this->checkNormalDomination(el, *succSuccPtr, *succPtr))
		{
			this->tableSequences.erase(succPtr->first);
			succPtr = this->LS.erase(succPtr);
//			std::cout << "Rimozione elemento DX skyline!\n";
		}
		else
			break;
	}
	// Verify whether "el" represents a stricter upper bound in the second dimension than the
	// current upper bound in the skyline: if so, make "el" the new upper bound.
	if(this->checkDomination(el, *succPtr))
	{
		this->tableSequences.erase(succPtr->first);
		this->LS.erase(succPtr);
	}

	// Return true (insertion successful).
	return(true);
}

/**
 * @brief This method merges "this" skyline with "LSOther".
 */
void LinearSkyline::mergeSkyline(LinearSkyline& LSOther)
{
	const auto ls = LSOther.getLS();
	for(const auto el : ls) this->insertElement(el, LSOther.getSequence(el.first));
}

void LinearSkyline::testLinSky()
{
	// Inizializza set punti random.
	std::srand(std::time(nullptr)); // use current time as seed for random generator
	// std::multiset<LinearSkyline::ElementSkyline> setPoints;
	std::vector<LinearSkyline::ElementSkyline> setPoints;
//	for(uint32_t i = 0; i < 500; i++)
//		setPoints.push_back(std::make_pair(std::rand() % 100, std::rand() % 300));
	setPoints.push_back(std::make_pair(27, 100));
	setPoints.push_back(std::make_pair(13, 158));
	setPoints.push_back(std::make_pair(32, 160));
	setPoints.push_back(std::make_pair(31, 98));
	setPoints.push_back(std::make_pair(26, 140));
	setPoints.push_back(std::make_pair(9, 160));
	setPoints.push_back(std::make_pair(22, 142));
	setPoints.push_back(std::make_pair(36, 158));

	// Trova gli upper-bound.
	auto it = setPoints.begin();
	LinearSkyline::ElementSkyline ub1 = *it; LinearSkyline::ElementSkyline ub2 = ub1; it++;
	for(; it != setPoints.end(); it++)
	{
		if((it->first < ub1.first) || (it->first == ub1.first && it->second < ub1.second)) ub1 = *it;
		if((it->second < ub2.second) || (it->second == ub2.second && it->first < ub2.first)) ub2 = *it;
	}
	std::cout << "Primo upper bound: " << ub1.first << "," << ub1.second << "\n";
	std::cout << "Secondo upper bound: " << ub2.first << "," << ub2.second << "\n";



	// Stampa elementi.
	for(auto el : setPoints)
		std::cout << "Val: " << el.first << "," << el.second << "\n";


	// 1 - Create and initialize a linear skyline with two upper bounds.
	LinearSkyline lin(ub1, std::vector<uint32_t>(0), ub2, std::vector<uint32_t>(0));

	// 2 - Remove the upper bounds from the ordered set.
	// setPoints.erase(ub1); setPoints.erase(ub2);

	// 3 - Filter out those points that lie in the area dominated by the upper bounds.
	// std::multiset<LinearSkyline::ElementSkyline> setFilteredPoints;
	std::vector<LinearSkyline::ElementSkyline> setFilteredPoints;
//	auto sub1 = lin.getFirstUpperBound(); // lin.cbegin();
//	auto sub2 = lin.getSecondUpperBound(); // ++lin.cbegin();
	for(it = setPoints.begin(); it != setPoints.end(); it++)
	{
//		if(*it == sub1 || *it == sub2) continue;
//		if(!(LinearSkyline::checkDomination(sub1, *it) || LinearSkyline::checkDomination(sub2, *it)))
			// setFilteredPoints.insert(*it);
			setFilteredPoints.push_back(*it);
	}

	// 4 - Check if the remaining points can be inserted into the linear skyline.
	std::cout << "Numero punti da valutare per inserzione skyline: " << setFilteredPoints.size() << "\n";
	for(auto el : setFilteredPoints)
	{
		std::cout << "Attempting insertion of (" << el.first << "," << el.second << "): ";
		std::cout << lin.insertElement(el) << "\n";
	}




	// Scrittura file test.
	std::ofstream myfile;
	myfile.open("points.txt");
	for(auto el : setPoints)
		myfile << el.first << " " << el.second << "\n";
	myfile.close();

	myfile.open("filteredPoints.txt");
	for(auto el : setFilteredPoints)
		myfile << el.first << " " << el.second << "\n";
	myfile.close();

	myfile.open("skyline.txt");
	for(auto it = lin.cbegin(); it != lin.cend(); it++)
		myfile << it->first << " " << it->second << "\n";
	// myfile << ub1.first << " " << ub1.second << "\n";
	// myfile << ub2.first << " " << ub2.second << "\n";
	myfile.close();
}
