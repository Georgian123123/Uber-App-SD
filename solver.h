// Copyright 2019 SD_Homework_Team
#ifndef SOLVER_H_
#define SOLVER_H_
#include <iostream>
#include <fstream>
#include <list>
#include <string>
#include <vector>
#include <queue>
#include <iomanip>
	// Hashtable.
	template<typename Tkey, typename Tvalue>
	// key & val
	struct elem_info {
		Tkey name;
		Tvalue val;
	};
	// hash function
	int myhash(std :: string str) {
		return str.length();
	}
	// Implemented Hashtable with class.
	template<typename Tkey, typename Tvalue>
	class Hashtable {
	private :
		std :: list<struct elem_info<Tkey, Tvalue> >*H;
	    int HMAX;
	    int size;
	    int (*hash)(Tkey);

	public :
	    	Hashtable() {}
			Hashtable(int &hmax, int(*h)(Tkey)) {
				HMAX = hmax;
				size = 0;
				this->hash = h;
				H = new std :: list<struct elem_info<Tkey, Tvalue> >[HMAX];
			}
			~Hashtable() {
				for (int i = 0; i < HMAX; i++){
					while (!H[i].empty())
						H[i].pop_front();
				}
				delete[] H;
			}
			// getting a value with a specific key.
	        void insert(Tkey key, Tvalue value) {
	           struct elem_info<Tkey, Tvalue> info;
	           int ok = 0;
	           int hkey = hash(key) % HMAX;
	           	for (auto it = this->H[hkey].begin();
	           		it != this->H[hkey].end(); ++it) {
	           		if(it->name == key) {
	           			it->val = value;
	           			ok = 1;
	           		}
	           	}
	           	if(ok == 0) {
	           		info.name = key;
	           		info.val = value;
	           		this->H[hkey].push_back(info);
	           		++size;
	           	}
	       }
	       // getting a key with a specific value.
	        Tvalue getValue(Tkey key) {
	            int hkey = hash(key) % HMAX;
	            for (auto it = this->H[hkey].begin();
	            	it != this->H[hkey].end(); ++it) {
	            	if(it->name == key) {
	            		return it->val;
	            	}
	            }
	           return -1;
			}
		};
		// Total uber players.
class Uber {
	private :
		int rating, totdist, totrides;
		std :: string nickname, station, state;
	public :
		Uber() {
			state = "off";
			rating = 0;
			totdist = 0;
			totrides = 0;
			nickname ="";
			station = "";
			state = "";
		}
		~Uber() {}
		void setNickname(std :: string str) {
			nickname = str;
		}
		void setStation(std :: string str) {
			station = str;
		}
		void setState(std :: string str) {
			state = str;
		}
		void setRating(int rate) {
			rating += rate;
		}
		std :: string getNickname() {
			return nickname;
		}
		std :: string getState() {
			return state;
		}
		std :: string getStation() {
			return station;
		}
		int getRating() {
			return rating;
		}
		float getMediumRat() {
			if(totrides == 0) {
				return 0;
			} else {
				return rating * 1.0 / totrides;
			}
		}
		void setDist(int x) {
			totdist += x;
		}
		int getDist() {
			return totdist;
		}
		void setRides(int x) {
			totrides += x;
		}
		int getRides() {
			return totrides;
		}
};
// our map
template <typename Tkey, typename Tvalue>
class Graph {
	private :
	std :: vector<std :: list<int>> node;
	int size;

	public :
	Hashtable<std :: string, int> *Map;
	std :: vector<std :: string> cities;
	explicit Graph(int n) {
		size = n;
		Map = new Hashtable<Tkey, Tvalue>(n, &myhash);
		for (int i = 0; i < n; i++) {
			std::list <int> list;
			node.push_back(list);
		}
	}
	Graph() {}
	~Graph() {
		delete Map;
	}
	// inserting a city
	void insert(Tkey key, Tvalue value) {
		Map -> insert(key, value);
	}
	// inserting a path between 2 cities
	void addEdge(Tkey departure, Tkey arrival) {
		int src = Map->getValue(departure);
		int dst = Map->getValue(arrival);
		if(src == dst) {
			return;
		}
		node[src].push_back(dst);
	}
	// removing a path between 2 cities
	void removeEdge(Tkey departure, Tkey arrival) {
		int src = Map->getValue(departure);
		int dst = Map->getValue(arrival);
		if(src == dst) {
			return;
		}
		for (auto it = node[src].begin(); it
			!= node[src].end(); it++) {
			if ((*it) == dst) {
				node[src].erase(it);
				return;
			}
		}
	}
	// verifing if it exists a direct path
	bool hasEdge(Tkey departure, Tkey arrival)  {
		int src = Map->getValue(departure);
		int dst = Map->getValue(arrival);
		if(src == dst) {
			return true;
		}
		for (auto it = node[src].begin();
			it != node[src].end(); it++) {
			if ((*it) == dst) {
				return true;
			}
		}
		return false;
	}
	// verifing if it exists a path
	bool hasPath(Tkey departure, Tkey arrival) const {
		int src = Map->getValue(departure);
		int dst = Map->getValue(arrival);
		if (src == dst) {
			return true;
		}
		std::queue <int> findPath;
		std::vector<bool> visited(size, 0);
		findPath.push(src);
		visited[src] = 1;
		while (!findPath.empty()) {
			int current = findPath.front();
			findPath.pop();
			for (auto it = node[current].begin();
				it != node[current].end(); it++){
				if (visited[*it] == 0) {
					findPath.push(*it);
					visited[*it] = 1;
				}
				if (*it == dst) {
					return true;
				}
			}
		}
		return false;
		}
		// getting all the neighbours
	std :: list<std :: string > hasPathforClient(
		std :: string& end, std :: list<std :: string> &list) {
		int src = Map->getValue(end);
		int i = 0;
		for (auto it = node[src].begin();
			it != node[src].end(); it++){
			i = 0;
			for(auto &it2 : cities) {
				if(*it == i) {
					list.push_back(it2);
				}
				i++;
			}
		}
		return list;
	}
	// getting costs for a path
    int hasPathlength(Tkey departure, Tkey arrival) const {
		int src = Map->getValue(departure);
		int dst = Map->getValue(arrival);
		if(src == dst) {
			return 0;
		}
		std::vector<bool> visited(size, 0);
		std::queue <int> findPath;
		int distances[size] = {0};
		distances[src] = 0;
		findPath.push(src);
		visited[src] = true;
		while(!findPath.empty()) {
			int current = findPath.front();
			findPath.pop();
			for (auto it = node[current].begin();
				it != node[current].end(); it++){
				if (visited[*it] == 0) {
					distances[*it] = distances[current] + 1;
					visited[*it] = 1;
					findPath.push(*it);
				}
			}
		}
		return distances[dst];
	}
};
class solver {
	private :
		Graph<std :: string, int> *graph;
		std :: vector<Uber> Ubers;
		int N, M;

	public:
	~solver() {
		delete graph;
	}
	void  task1_solver(std::ifstream& fin, std::ofstream& fout) {
		// added the cities and the paths between them
		// verifing if it exists the path
		std :: string city;
		fin >> N >> M;
		graph = new Graph<std :: string, int>(N);
		for(int i = 0; i < N; ++i) {
			fin >> city;
			graph->insert(city, i);
			graph->cities.push_back(city);
		}
		std :: string city1, city2;
		for(int i = 0; i < M; ++i) {
			fin >> city1 >> city2;
			graph->addEdge(city1, city2);
		}
		int querry;
		fin >> querry;
		for(int i = 0; i < querry; ++i) {
			fin >> city1 >> city2;
			if(graph->hasPath(city1, city2)) {
				fout << "y" << std :: endl;
			} else {
				fout << "n" << std :: endl;
			}
		}
	}
	void task2_solver(std::ifstream& fin, std::ofstream& fout) {
		// verify if it exists the path & show the cost
		int querry;
		fin >> querry;
		std :: string city1, city2;
		for(int i = 0; i < querry; ++i) {
			fin >> city1 >> city2;
			if(graph->hasPath(city1, city2)) {
				fout << graph->hasPathlength(city1, city2) << std :: endl;
			} else {
				fout << -1 << std :: endl;
			}
		}
	 }

		void task3_solver(std::ifstream& fin, std::ofstream& fout) {
			// adding the paths, removing them, and make the path bidirectional
			std :: string TODO, city1, city2;
			int querry, number;
			fin >> querry;
			for (int i = 0; i < querry; ++i) {
				fin >> TODO;
				if(TODO == "c") {
					fin >> city1 >> city2;
					fin >> number;
					switch(number) {
						case 0 :
								if(!graph->hasEdge(city1, city2)) {
									graph->addEdge(city1, city2);
								}
							break;
						case 1 :
								if(graph->hasEdge(city2, city1)) {
									graph->removeEdge(city2, city1);
								}
								if(graph->hasEdge(city1, city2)) {
								graph->removeEdge(city1, city2);
								}
								break;
						case 2 :
							if (!graph->hasEdge(city1, city2)) {
								graph->addEdge(city1, city2);
							}
							if (!graph->hasEdge(city2, city1)) {
								graph->addEdge(city2, city1);
							}
							break;
						case 3 :
							if (graph->hasEdge(city1, city2) == 1
								&& graph->hasEdge(city2, city1) == 0) {
								graph->removeEdge(city1, city2);
								graph->addEdge(city2, city1);
							} else if (graph->hasEdge(city1, city2) == 0
								&& graph->hasEdge(city2, city1) == 1) {
								graph->removeEdge(city2, city1);
								graph->addEdge(city1, city2);
							}
							break;
					}
				}
				// show if it exists the paths.
				if(TODO == "q") {
					fin >> city1 >> city2;
					fin >> number;
					switch(number) {
						case 0 :
							if(graph->hasPath(city1, city2)) {
								fout << "y" << std :: endl;
							} else {
								fout << "n" << std :: endl;
							}
							break;
						case 1 :
							if(graph->hasPath(city1, city2)) {
								fout << graph->hasPathlength(city1, city2)
								<< std :: endl;
							} else {
								fout << -1 << std :: endl;
							}
							break;
						case 2 :
							std :: string city3;
							fin >> city3;
							  if (graph->hasPath(city1, city3) &&
							  	graph->hasPath(city3, city2)) {
									fout <<
								(graph->hasPathlength(city1, city3) +
									graph->hasPathlength(city3, city2))
									<< std :: endl;
								} else {
									fout << -1 << std :: endl;
								}
							break;
					}
				}
			}
		}
void sortingBubble(float x[], std :: string y[], int &max) {
	// sorting function for displays alike : top (task4)
    for (int i = 0; i < max - 1; ++i) {
    	for(int j = i + 1; j < max; ++j) {
    		if(x[i] < x[j]) {
    			float aux = x[i];
    			x[i] = x[j];
    			x[j] = aux;
    			std :: string aux2 = y[i];
	    		y[i] = y[j];
	    		y[j] = aux2;
    		} else if (x[i] == x[j]) {
    					if (y[i] > y[j]) {
    						float aux3 = x[i];
    						x[i] = x[j];
    						x[j] = aux3;
	    					std :: string aux2 = y[i];
	    					y[i] = y[j];
	    					y[j] = aux2;
	    				}
    				}
    			}
    		}
	}
void getUber(int distances[], std :: string names[], float rats[], int max) {
	// getting the nearest driver comparing the distances, rating, and names
	// for all
	for(int i = 0; i < max - 1; ++i) {
		for(int j =  i + 1; j < max ; ++j) {
			if(distances[i] > distances[j]) {
				int aux = distances[i];
				distances[i] = distances[j];
				distances[j] = aux;
				std :: string aux2 = names[i];
				names[i] = names[j];
				names[j] = aux;
				float x = rats[i];
				rats[i] = rats[j];
				rats[j] = x;
			} else {
				if(distances[i] == distances[j]) {
					if(rats[i] < rats[j]) {
						int aux3 = distances[i];
						distances[i] = distances[j];
						distances[j] = aux3;
						std :: string aux4 = names[i];
						names[i] = names[j];
						names[j] = aux4;
						float y = rats[i];
						rats[i] = rats[j];
						rats[j] = y;
					} else {
						if(rats[i] == rats[j]) {
							if(names[i] > names[j]) {
								int aux5 = distances[i];
								distances[i] = distances[j];
								distances[j] = aux5;
								std :: string aux6 = names[i];
								names[i] = names[j];
								names[j] = aux6;
								float z = rats[i];
								rats[i] = rats[j];
								rats[j] = z;
							}
						}
					}
				}
			}
		}
	}
}
void task4_solver(std::ifstream& fin, std::ofstream& fout) {
	int querry;
	fin >> querry;
	std :: string TODO;
	std :: string name, city;
	for(int m = 0; m < querry; ++m) {
		// adding the drivers
		fin >> TODO;
		if (TODO == "d") {
			int ok = 0;
			fin >> name >> city;
			if(Ubers.size() == 0) {
				Uber player1;
				player1.setState("online");
				player1.setStation(city);
				player1.setNickname(name);
				Ubers.push_back(player1);
			} else {
				for(auto &it : Ubers) {
					if(it.getNickname() == name) {
						if(it.getState() == "off") {
							it.setState("online");
							it.setStation(city);
							ok = 1;
							break;
						} else {
							it.setStation(city);
							ok = 1;
							break;
						}
					}
				}
					if(ok == 0) {
						Uber player;
						player.setState("online");
						player.setStation(city);
						player.setNickname(name);
						Ubers.push_back(player);
					}
				}
		} else {
		if (TODO == "b") {
			// make driver offline
			std :: string name;
			fin >> name;
			for(auto &it : Ubers) {
				if(it.getNickname() == name) {
					it.setState("off");
						break;
				}
			}
		} else {
		if (TODO == "r") {
			std :: string start, end;
			int rating , ok = 0;
			fin >> start >> end >> rating;
			Uber minUber;
			int distances[Ubers.size()];
			float rats[Ubers.size()];
			std :: string names[Ubers.size()];
			int i = 0;
		for(auto &it : Ubers) {
			// find if it exists a driver able to reach the client
			if(graph->hasPath(it.getStation(), start)
				&& it.getState() == "online") {
				distances[i] = graph->hasPathlength(it.getStation(), start);
				names[i] = it.getNickname();
				rats[i] = it.getMediumRat();
				ok = 1;
				++i;
			}
		}
		if(ok == 0) {
			// it doesn't
			fout << "Soferi indisponibili" << std :: endl;
		} else {
			// got the driver.Verify if it's possible to arrive at the
			// destination.If not, arrive to the first neighbour.
			getUber(distances, names, rats,  i);
			for(auto &it : Ubers) {
				if(it.getNickname() == names[0]) {
					minUber = it;
					break;
				}
			}
			int op = 0;
			int x = graph->hasPathlength(minUber.getStation(), start);
			for(auto &it : Ubers) {
				if(it.getNickname() == minUber.getNickname()) {
					if(graph->hasPath(start, end)) {
						it.setStation(end);
						it.setRating(rating);
						it.setRides(1);
						it.setDist(graph->hasPathlength(start, end) + x);
						op = 1;
						break;
					} else {
						std :: list <std :: string> list;
						list = graph->hasPathforClient(end, list);
						if(list.size() == 0) {
							// it doesn't exist a neightbour so display it.
							fout << "Destinatie inaccesibila" << std :: endl;
							op = 1;
							break;
						}
					   	for (auto &it2 : list) {
						if(graph->hasPath(start, it2)) {
							it.setStation(it2);
							it.setRating(rating);
							it.setRides(1);
							it.setDist(graph->hasPathlength(start, it2) + x);
							op = 1;
							break;
						}
					  }
				   }
				}
			}
			if(op == 0) {
				// it doesn't exist a neightbour so display it.
				fout << "Destinatie inaccesibila" << std :: endl;
			}
		}
	} else {
		if(TODO == "top_rating") {
			// display the requested number of drivers.
			// they are sorted by rating and name.
			size_t numb;
			fin >> numb;
			if(Ubers.size() != 0) {
				float rat_m2[Ubers.size()];
				std :: string names[Ubers.size()];
				int i = 0;
				for(auto &it : Ubers) {
					rat_m2[i] = it.getMediumRat();
					names[i] = it.getNickname();
					++i;
				}
				size_t min;
				if(numb < Ubers.size()) {
					min = numb;
				} else {
					min = Ubers.size();
				}
				sortingBubble(rat_m2, names, i);
				for(size_t i = 0; i < min; ++i) {
					fout << names[i]<< ":" <<std :: fixed <<
					std :: setprecision(3) << rat_m2[i] << " ";
				}
				fout << std :: endl;
			} else {
				fout << std :: endl;
			}
		} else {
		if (TODO == "top_dist")	{
			// display the requested number of drivers.
			// they are sorted by total distance and name.
			size_t numb;
			fin >> numb;
			if(Ubers.size() != 0) {
				float tot_dist[Ubers.size()];
				std :: string names[Ubers.size()];
				int i = 0;
				for(auto &it : Ubers) {
					tot_dist[i] = it.getDist();
					names[i] = it.getNickname();
					++i;
				}
				size_t min;
				if(numb < Ubers.size()) {
					min = numb;
				} else {
					min = Ubers.size();
				}
				sortingBubble(tot_dist, names, i);
				for(size_t i = 0; i < min; ++i) {
					fout << names[i]<< ":" <<
					std :: fixed << std :: setprecision(0) << tot_dist[i] << " ";
				}
				fout << std :: endl;
			} else {
				fout << std :: endl;
			}
		} else {
		if(TODO == "top_rides") {
			// display the requested number of drivers.
			// they are sorted by total rides and name.
			size_t numb;
			fin >> numb;
			float tot_rides[Ubers.size()];
			std :: string names[Ubers.size()];
			if(Ubers.size() != 0) {
				int i = 0;
				for(auto &it : Ubers) {
					tot_rides[i] = it.getRides();
					names[i] = it.getNickname();
					++i;
				}
				size_t min;
				if(numb < Ubers.size()) {
					min = numb;
				} else {
					min = Ubers.size();
				}
				sortingBubble(tot_rides, names, i);
				for(size_t i = 0; i < min; ++i) {
					fout << names[i]<< ":" <<
					std :: fixed << std :: setprecision(0)<< tot_rides[i] << " ";
				}
				fout << std :: endl;
			} else {
				fout << std :: endl;
			}
		} else {
			if (TODO == "info") {
				// display all the info about a driver
				std :: string name;
				fin >> name;
				if(Ubers.size() != 0) {
					for (auto &it : Ubers) {
						if(it.getNickname() == name) {
								float x = it.getMediumRat();
								fout << it.getNickname() << ": " <<
								it.getStation() << " " <<std :: fixed <<
								std :: setprecision(3) << x  << " " <<
								it.getRides() << " " <<it.getDist()
								<< " " << it.getState() << std :: endl;
						}
					}
				} else {
					fout << std :: endl;
				}
			}
		}
						}
					}
				}
			}
		}
	}
}
void task5_solver(std::ifstream& fin, std::ofstream& fout) {
	int dist, nr_cities;
	std :: string name, city;
	fin >> dist >> name >> nr_cities;
	std :: vector<std :: string > cities;
	std :: vector<std :: string > cities2;
	std :: vector<int> dists;
	for(int i = 0; i < nr_cities; ++i) {
		fin >> city;
		cities.push_back(city);
	}
	int ok = 0;
	// verify if it exists path and the cost.
	// verify if the gas it's enough
	// if yes, put it in a vector and after that, sorted them.
	for(auto &it : Ubers) {
		if(it.getNickname() == name) {
			for(auto &it2 : cities) {
				ok = 0;
				if(graph->hasPath(it.getStation(), it2)) {
					if(dist >= graph->hasPathlength(it.getStation(),
						it2)) {
						for(auto &it3 : cities2) {
							if(it3 == it2) {
								ok = 1;
							}
						}
						if(ok == 0) {
							cities2.push_back(it2);
							dists.push_back(
								graph->hasPathlength(it.getStation(),
								it2));
						}
					}
				}
			}
		}
	}		int kDim = cities2.size();
			std :: string names[kDim];
			int distances[kDim];
			int i = 0;
			for(auto &it : cities2) {
				names[i] = it;
				i++;
			}
			i = 0;
			for(auto &it : dists) {
				distances[i] = it;
				i++;
			}
			for(int j = 0; j < i - 1; ++j) {
				for(int k = j + 1; k < i; ++k) {
					if(distances[j] > distances[k]) {
						int aux = distances[j];
						distances[j] = distances[k];
						distances[k] = aux;
						std :: string aux3;
						aux3 = names[j];
						names[j] = names[k];
						names[k] = aux3;
					} else if (distances[j] == distances[k]) {
						if(names[j] > names[k]) {
							int aux4 = distances[j];
							distances[j] = distances[k];
							distances[k] = aux4;
							std :: string aux2 = names[j];
							names[j] = names[k];
							names[k] = aux2;
						}
					}
				}
			}
			for(int j = 0; j < i; ++j) {
				fout << names[j] << " ";
			}
		}
};
#endif  // SOLVER_H_
