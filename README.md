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

Please use the ```help``` option when executing the application binary to see all the possible parameters that you can pass to the application. The application have two different modalities:

- Compute TASeR queries: this requires to pass several parameters, i.e.,
  - the files of the vertices and edges of a road network (see the -v and -e parameters);
  - the POIs of COIs of choice (see the -p parameters);
  - the approach used to compute the queries (see the -a parameter)
  - file containing the precoumpted distances between POIs of COIs to be used to speed up the queries' computation (see the -l parameter)
  - the length of the COI sequences of the TASeR queries (see the -m parameter)
  - the number of repetitions (see the -r parameter)
    
- Precompute distances between POIs of COIs: in this mode, which is triggered when using the `````` parameter, the application takes as input the vertices and edges of a road network, plus the POIs of COIs of choice, and precompute the distances between the POIs over the road netork.

Below is an example of how to execute our application to compute TASeR queries. 

```
./bin/CASRQ -v datasets/Oslo/roadnetwork/RoadVerticesOSLO.txt
            -e datasets/Oslo/roadnetwork/RoadEdgesOSLO.txt 
            -p datasets/Oslo/poi/uniform/AtmBanksOSLO.txt.CORRECTED.txt.UNIFORM.txt
            -p datasets/Oslo/poi/uniform/CoffeeShopsOSLO.txt.CORRECTED.txt.UNIFORM.txt
            -p datasets/Oslo/poi/uniform/GasStationsOSLO.txt.CORRECTED.txt.UNIFORM.txt
            -p datasets/Oslo/poi/uniform/MovieTheatersOSLO.txt.CORRECTED.txt.UNIFORM.txt
            -p datasets/Oslo/poi/uniform/PharmaciesOSLO.txt.CORRECTED.txt.UNIFORM.txt
            -p datasets/Oslo/poi/uniform/PubsBarsOSLO.txt.CORRECTED.txt.UNIFORM.txt
            -p datasets/Oslo/poi/uniform/RestaurantsOSLO.txt.CORRECTED.txt.UNIFORM.txt
            -a 0 -m 2 -r 10
```


## Cite us

Please cite our MDM 2020 article if you have found our contributions useful, or you have used them within your work.

```bibtex
@inproceedings{lettich2020,
  author={Lettich, Francesco and Nascimento, Mario and Anwar, Samiul},
  booktitle={2020 21st IEEE International Conference on Mobile Data Management (MDM)}, 
  title={Trade-off Aware Sequenced Routing Queries (or OSR Queries when POIs are not Free)}, 
  year={2020},
  volume={},
  number={},
  pages={59-68},
  keywords={Roads;Routing;Task analysis;Upper bound;Conferences;Indexes;Fuels;Optimal Sequenced Routing;Road Networks;Linear Skylines;Trade-off Aware Sequenced Routing},
  doi={10.1109/MDM48529.2020.00027}}
