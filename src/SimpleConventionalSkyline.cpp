// *** INCLUDES *** //

#include "SimpleConventionalSkyline.hpp"



// *** PUBLIC CTORS DEFINITIONS *** //

/**
 * @brief Initializes a new linear skyline by inserting the lower bound associated
 *  	  with the first type of cost.
 *
 * @note We are assuming that we know the two upper bounds and elements are inserted into
 * 	     the skyline in increasing order of the first type of cost.
 */
SimpleConventionalSkyline::SimpleConventionalSkyline(const ElementSkyline& ub1, const std::vector<uint32_t>& ub1Seq) :
SimpleSkyline(ub1, ub1Seq)
{}



// *** PUBLIC METHODS DEFINITIONS *** //

/**
 * @brief This method attempts to insert "el" into the linear skyline, along with the associated sequence
 * 		  of POIs.
 *
 * @return true if the insertion is successful, false otherwise.
 * @note The method may remove elements from the linear skyline if "el" is inserted into the skyline.
 */
bool SimpleConventionalSkyline::insertElement(const ElementSkyline& el, const std::vector<uint32_t>& elSeq)
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


	// Insert el and return true.
	this->LS.push_back(el);
	this->LSS.push_back(elSeq);
	return(true);
}
