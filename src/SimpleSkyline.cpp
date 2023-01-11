/**
 * @brief Class modeling a linear skyline.
 * @author Francesco Lettich
 */

// *** INCLUDES *** //

#include "SimpleSkyline.hpp"

#include <vector>
#include <set>
#include <ctime>
#include <iostream>
#include <cstdlib>
#include <fstream>



// *** PUBLIC CTORS DEFINITIONS *** //

SimpleSkyline::SimpleSkyline(const ElementSkyline& ub1)
{
	this->LS.push_back(ub1);
}

/**
 * @brief Initializes a new linear skyline by inserting the lower bound associated
 *  	  with the first type of cost.
 *
 * @note We are assuming that we know the two upper bounds and elements are inserted into
 * 	     the skyline in increasing order of the first type of cost.
 */
SimpleSkyline::SimpleSkyline(const ElementSkyline& ub1, const std::vector<uint32_t>& ub1Seq)
{
	this->LS.push_back(ub1);
	this->LSS.push_back(ub1Seq);
}



// *** PUBLIC METHODS DEFINITIONS *** //

/**
 * @brief This method attempts to insert "el" into the linear skyline, along with the associated sequence
 * 		  of POIs.
 *
 * @return true if the insertion is successful, false otherwise.
 * @note The method may remove elements from the linear skyline if "el" is inserted into the skyline.
 */
bool SimpleSkyline::insertElement(const ElementSkyline& el, const std::vector<uint32_t>& elSeq)
{
	uint32_t indexLeftNeigh = this->LS.size() - 1;
//	std::cout << "Left neighbor: " << *lneigh.first << "," << lneigh.second << "\n";


	// First, check if the left neighbor conventionally dominates "el": if so, refuse the insertion.
	// NOTE: we discard also those routes that have travel distance and travel cost equal to
	// 		 the left neighbor. This is compatible with the definition of linear skylines being
	//		 *ordered* sets.
	// if(this->checkDomination(this->LS[indexLeftNeigh], el)) return(false);
	if(this->LS[indexLeftNeigh].second <= el.second) return(false);
	//	std::cout << "Non convenzionalmente dominato.\n";

	// Manage the case where the new element has the same travel cost but smaller travel distance
	// than its left neighbor: remove the left neighbor from the skyline and update the index accordingly.
	if(this->LS[indexLeftNeigh].first == el.first)
	{
		this->LS.pop_back();
		this->LSS.pop_back();
		indexLeftNeigh -= (indexLeftNeigh > 0);
	}


	// "el" can be inserted. Verify if there are skyline elements on the left of "el"
	// that are linearly dominated.
//	std::cout << "Verifica se ci sono elementi da rimuovere a SX nella skyline.\n";
	while(indexLeftNeigh > 0)
	{
		if(this->checkNormalDomination(this->LS[indexLeftNeigh - 1], el, this->LS[indexLeftNeigh]))
		{
			this->LS.pop_back();
			this->LSS.pop_back();
//			std::cout << "Rimozione elemento SX skyline!\n";
		}
		else
			break;

		indexLeftNeigh--;
	}


	// Insert el and return true.
	this->LS.push_back(el);
	this->LSS.push_back(elSeq);
	return(true);
}

/**
 * @brief Internal test method to verify that this class works correctly.
 */
void SimpleSkyline::testLinSky()
{
	// Inizializza set punti random.
	std::srand(std::time(nullptr)); // use current time as seed for random generator
	std::multiset<SimpleSkyline::ElementSkyline> setPoints;
	// std::vector<SimpleSkyline::ElementSkyline> setPoints;
//	for(uint32_t i = 0; i < 500; i++)
//		setPoints.insert(std::make_pair(std::rand() % 500, std::rand() % 500));
	setPoints.insert(std::make_pair(27, 100));
	setPoints.insert(std::make_pair(13, 158));
	setPoints.insert(std::make_pair(32, 160));
	setPoints.insert(std::make_pair(31, 98));
	setPoints.insert(std::make_pair(26, 140));
	setPoints.insert(std::make_pair(9, 160));
	setPoints.insert(std::make_pair(22, 142));
	setPoints.insert(std::make_pair(36, 158));

	// Trova gli upper-bound.
	auto it = setPoints.begin();
	SimpleSkyline::ElementSkyline ub1 = *it; SimpleSkyline::ElementSkyline ub2 = ub1; it++;
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


	// 1 - Create and initialize a linear skyline.
	SimpleSkyline lin(ub1);

	// 2 - Remove the upper bounds from the ordered set.
	setPoints.erase(ub1); // setPoints.erase(ub2);

	// 3 - Filter out those points that lie in the area dominated by the upper bounds.
	// std::multiset<LinearSkyline::ElementSkyline> setFilteredPoints;
	std::vector<SimpleSkyline::ElementSkyline> setFilteredPoints;
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

/**
 * @brief This method merges this linear skyline with another one passed via reference.
 */
SimpleSkyline SimpleSkyline::mergeWithSkyline(const SimpleSkyline& SS)
{
	// TODO: da testare.
	auto it = this->LS.cbegin(); auto itSeq = this->LSS.cbegin();
	auto itOther = SS.LS.cbegin(); auto itOtherSeq = SS.LSS.cbegin();
	SimpleSkyline res(*it, *itSeq); // Pull in the upper bound in the first dimension.

	// Perform a 2-way linear scan of the two linear skylines.
	it++; itSeq++; itOther++; itOtherSeq++;
	while(it != this->LS.cend() && itOther != SS.LS.cend())
	{
		if(it->first <= it->second)
		{
			res.insertElement(*it, *itSeq);
			it++; itSeq++;
		}
		else
		{
			res.insertElement(*itOther, *itOtherSeq);
			itOther++; itOtherSeq++;
		}
	}

	// Add remaining elements in one of the two skylines.
	while(it != this->LS.cend())
	{
		res.insertElement(*it, *itSeq);
		it++; itSeq++;
	}
	while(itOther != SS.LS.cend())
	{
		res.insertElement(*itOther, *itOtherSeq);
		itOther++; itOtherSeq++;
	}

	return(res);
}
