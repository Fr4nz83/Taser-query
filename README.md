# Trade-off aware sequenced routing queries (or OSR queries when pois are not free)

This repository contains the code behind the paper *Trade-off aware sequenced routing queries (or OSR queries when pois are not free)*, by Francesco Lettich, Mario A. Nascimento, and Samiul Anwar; all the authors were working at the University of Alberta when working on this paper.

The paper has been accepted and presented at the [21st IEEE International Conference on Mobile Data Management (MDM 2020)](https://mdmconferences.org/mdm2020/).


## How to compile the code

The code behind the paper has been written in C++. You are expected to compile it under some linux distribution, e.g., Ubuntu. 
In order to correctly compile and execute the application, the following shared libraries need to be present (they can be installed via ```apt install```):

- Boost
- [Intel TBB](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onetbb.html)

  We have also included the following libraries under the ```lib``` folder, which are also required to compile the code:
  
- [GeographicLib](https://geographiclib.sourceforge.io/)
- [Lemon Graph Library](https://en.wikipedia.org/wiki/LEMON_(C%2B%2B_library))

To compile the code, simply execute ```make all``` within the project's root folder.


## Execution example

Below is an example of how to execute our application. Please ... 

./bin/CASRQ -v datasets/Oslo/roadnetwork/RoadVerticesOSLO.txt -e datasets/Oslo/roadnetwork/RoadEdgesOSLO.txt -p datasets/Oslo/poi/uniform/AtmBanksOSLO.txt.CORRECTED.txt.UNIFORM.txt -p datasets/Oslo/poi/uniform/CoffeeShopsOSLO.txt.CORRECTED.txt.UNIFORM.txt -p datasets/Oslo/poi/uniform/GasStationsOSLO.txt.CORRECTED.txt.UNIFORM.txt -p datasets/Oslo/poi/uniform/MovieTheatersOSLO.txt.CORRECTED.txt.UNIFORM.txt -p datasets/Oslo/poi/uniform/PharmaciesOSLO.txt.CORRECTED.txt.UNIFORM.txt -p datasets/Oslo/poi/uniform/PubsBarsOSLO.txt.CORRECTED.txt.UNIFORM.txt datasets/Oslo/poi/uniform/RestaurantsOSLO.txt.CORRECTED.txt.UNIFORM.txt -a 0 -m 2 -r 10


## Cite us


