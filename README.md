# Trade-off aware sequenced routing queries (or OSR queries when pois are not free)

This repository contains the code behind the paper *Trade-off aware sequenced routing queries (or OSR queries when pois are not free)* by Francesco Lettich, Mario A. Nascimento, and Samiul Anwar; the authors were all working at the University of Alberta when working at this paper.

The paper has been accepted and presented at the 21st IEEE International Conference on Mobile Data Management (MDM 2020).

## How to compile the code

The code has been written in C++. In order to correctly compile and execute the application, the following shared libraries need to be present:

- Boost
- Intel TBB
- GeographicLib

## Example of a single run 

./bin/CASRQ -v datasets/Oslo/roadnetwork/RoadVerticesOSLO.txt -e datasets/Oslo/roadnetwork/RoadEdgesOSLO.txt -p datasets/Oslo/poi/uniform/AtmBanksOSLO.txt.CORRECTED.txt.UNIFORM.txt -p datasets/Oslo/poi/uniform/CoffeeShopsOSLO.txt.CORRECTED.txt.UNIFORM.txt -p datasets/Oslo/poi/uniform/GasStationsOSLO.txt.CORRECTED.txt.UNIFORM.txt -p datasets/Oslo/poi/uniform/MovieTheatersOSLO.txt.CORRECTED.txt.UNIFORM.txt -p datasets/Oslo/poi/uniform/PharmaciesOSLO.txt.CORRECTED.txt.UNIFORM.txt -p datasets/Oslo/poi/uniform/PubsBarsOSLO.txt.CORRECTED.txt.UNIFORM.txt datasets/Oslo/poi/uniform/RestaurantsOSLO.txt.CORRECTED.txt.UNIFORM.txt -a 0 -m 2 -r 10
