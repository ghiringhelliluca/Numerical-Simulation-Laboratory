/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <math.h>
#include <set>
#include <vector>
#include <fstream>
#include "individual.h"

using namespace std;
using namespace arma;

void Individual::initialize(Random &_rnd, const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Cannot open file " << filename << endl;
        exit(EXIT_FAILURE);
    }
    
    vector<double> longitudes, latitudes;
    double lon, lat;
    
    // Read coordinates from file
    while (file >> lon >> lat) {
        longitudes.push_back(lon);
        latitudes.push_back(lat);
    }
    file.close();
    
    _ncity = longitudes.size();
    if (_ncity == 0) {
        cerr << "Error: No data found in file " << filename << endl;
        exit(EXIT_FAILURE);
    }
    
    // Initialize base individual with coordinates from file
    _base_individual.set_size(_ncity);
    _individual.set_size(_ncity);
    
    for (int i = 0; i < _ncity; i++) {
        _base_individual(i).resize(_ndim + 1);
        _base_individual(i)(0) = i + 1; // province ID
        _base_individual(i)(1) = longitudes[i]; // longitude
        _base_individual(i)(2) = latitudes[i];  // latitude
        
        _individual(i) = _base_individual(i); // copy to individual
    }
    

}

void Individual::initialize_from_base(Random &_rnd) {
    _ncity = _base_individual.n_elem;
    _individual.set_size(_ncity);
    
    for (int i = 0; i < _ncity; i++) {
        _individual(i) = _base_individual(i);
    }
}

// Haversine formula for calculating distance between two points on Earth
double Individual::haversine_distance(double lon1, double lat1, double lon2, double lat2) const {
    const double R = 6371.0; // Earth's radius in kilometers
    
    // Convert degrees to radians
    double lat1_rad = lat1 * M_PI / 180.0;
    double lat2_rad = lat2 * M_PI / 180.0;
    double dlat = (lat2 - lat1) * M_PI / 180.0;
    double dlon = (lon2 - lon1) * M_PI / 180.0;
    
    double a = sin(dlat/2) * sin(dlat/2) + 
               cos(lat1_rad) * cos(lat2_rad) * 
               sin(dlon/2) * sin(dlon/2);
    double c = 2 * atan2(sqrt(a), sqrt(1-a));
    
    return R * c;
}

// permutation function, make a pair permutation fix first province
void Individual::permutation(Random &_rnd) {
    for (int i = 1; i < _ncity; i++) {
        int j = int(_rnd.Rannyu(1, _ncity));
        if (i != j) {
            vec temp = _individual(i);
            _individual(i) = _individual(j);
            _individual(j) = temp;
        }
    }
}



double Individual::loss() {
    double total_distance = 0.0;
    
    for (int i = 0; i < _ncity; i++) {
        int next_i = (i == _ncity - 1) ? 0 : i + 1;
        
        double lon1 = _individual(i)(1);
        double lat1 = _individual(i)(2);
        double lon2 = _individual(next_i)(1);
        double lat2 = _individual(next_i)(2);
        
        total_distance += haversine_distance(lon1, lat1, lon2, lat2);
    }
    
    return total_distance;
}


pair<Individual, Individual> Individual::crossover_pair(const Individual& partner, Random& _rnd) {
    Individual mom = *this;
    Individual dad = partner;

    vector<int> order1, order2;
    int cut = int(_rnd.Rannyu(1, _ncity));

    // mother: this
    for (int i = 0; i < cut; ++i) {
        order1.push_back(mom._individual(i)(0));
    }

    for (int i = 0; i < _ncity; ++i) {
        int id = dad._individual(i)(0);
        if (find(order1.begin(), order1.end(), id) == order1.end()) {
            order1.push_back(id);
        }
    }

    // father: partner
    for (int i = 0; i < cut; ++i)
        order2.push_back(dad._individual(i)(0));

    for (int i = 0; i < _ncity; ++i) {
        int id = mom._individual(i)(0);
        if (find(order2.begin(), order2.end(), id) == order2.end())
            order2.push_back(id);
    }

    Individual child1, child2;
    child1._individual.set_size(_ncity);
    child2._individual.set_size(_ncity);
    child1._base_individual = _base_individual;
    child2._base_individual = _base_individual;
    child1._ncity = _ncity;
    child2._ncity = _ncity;

    // build child1._individual from order1
    for (int i = 0; i < _ncity; ++i) {
        int province_id = order1[i];
        for (int j = 0; j < _ncity; ++j) {
            if (_base_individual(j)(0) == province_id) {
                child1._individual(i) = _base_individual(j);
                break;
            }
        }
    }

    // build child2._individual from order2
    for (int i = 0; i < _ncity; ++i) {
        int province_id = order2[i];
        for (int j = 0; j < _ncity; ++j) {
            if (_base_individual(j)(0) == province_id) {
                child2._individual(i) = _base_individual(j);
                break;
            }
        }
    }

    return make_pair(child1, child2);
}

void Individual::mutation_shift(Random &_rnd) {
    int max_start = _ncity - 2;
    if (max_start < 1) return;

    int m = int(_rnd.Rannyu(1, max_start + 1));
    int max_block_len = min(4, _ncity - m);
    int block_len = int(_rnd.Rannyu(1, max_block_len + 1));

    int max_shift = _ncity - m - block_len;
    if (max_shift < 1) return;

    int shift = int(_rnd.Rannyu(1, max_shift + 1));

    vector<vec> new_path;

    // Before block
    for (int i = 0; i < m; ++i) {
        new_path.push_back(_individual(i));
    }
    // Between block and shift
    for (int i = m + block_len; i < m + block_len + shift; ++i) {
        new_path.push_back(_individual(i));
    }
    // Insert shifted block
    for (int i = 0; i < block_len; ++i) {
        new_path.push_back(_individual(m + i));
    }
    // Final part
    for (int i = m + block_len + shift; i < _ncity; ++i) {
        new_path.push_back(_individual(i));
    }

    for (int i = 0; i < _ncity; ++i) {
        _individual(i) = new_path[i];
    }
}

void Individual::mutation_inversion(Random &_rnd) {
    int max_start = _ncity - 3;
    int start = int(_rnd.Rannyu(1, max_start + 1));
    int end = int(_rnd.Rannyu(start + 1, _ncity));

    while (start < end) {
        vec tmp = _individual(start);
        _individual(start) = _individual(end);
        _individual(end) = tmp;
        ++start;
        --end;
    }
}

void Individual::mutation_swap_blocks(Random &_rnd) {
    int m = int(_rnd.Rannyu(1, int(_ncity / 2)));
    int a = int(_rnd.Rannyu(1, _ncity - 2 * m));
    int b = int(_rnd.Rannyu(a + m, _ncity - m + 1));

    for (int i = 0; i < m; ++i) {
        vec tmp = _individual(a + i);
        _individual(a + i) = _individual(b + i);
        _individual(b + i) = tmp;
    }
}

bool Individual::check() {
    // Check if first province is ID 1
    if (_individual(0)(0) != 1) {
        cout << "Error: first province is not ID 1!" << endl;
        return false;
    }

    // Check all provinces appear only once
    set<int> seen;
    for (int i = 0; i < _ncity; ++i) {
        int province_id = _individual(i)(0);
        if (seen.find(province_id) != seen.end()) {
            cout << "Error: province " << province_id << " appears more than once!" << endl;
            return false;
        }
        seen.insert(province_id);
    }

    // Check all provinces from 1 to n are present
    for (int i = 1; i <= _ncity; ++i) {
        if (seen.find(i) == seen.end()) {
            cout << "Error: missing province " << i << endl;
            return false;
        }
    }

    return true;
}

void Individual::print_path() {
    for (int i = 0; i < _ncity; ++i) {
        cout << _individual(i)(0) << " ";
    }
    cout << endl;
}

void Individual::save_path(ofstream& path) {
    for (int i = 0; i < _ncity; ++i) {
        path << _individual(i)(0) << " ";
    }
    path << endl;
}

void Individual::save_individual() {
    ofstream pos_individual("position_provinces.dat");
    pos_individual << "PROVINCE_ID      LONGITUDE      LATITUDE" << endl;

    for (int i = 0; i < _ncity; i++) {
        pos_individual << _individual(i)(0) << "             " 
                      << _individual(i)(1) << "             " 
                      << _individual(i)(2) << endl;
    }
    pos_individual.close();
}

vector<int> Individual::encode() const {
    vector<int> data;
    for (int i = 0; i < _ncity; ++i) {
        data.push_back(_individual(i)(0)); // province ID
    }
    return data;
}

void Individual::decode(const vector<int>& data) {
    if (data.empty()) return;
    
    _ncity = data.size();
    _individual.set_size(_ncity);
    
    for (int i = 0; i < _ncity; ++i) {
        _individual(i).resize(_ndim + 1);
        int province_id = data[i];
        _individual(i)(0) = province_id;
        _individual(i)(1) = GetXCoordinate(province_id);
        _individual(i)(2) = GetYCoordinate(province_id);
    }
}

// wrong function
double Individual::GetXCoordinate(int province_id) const {

    for (int i = 0; i < _base_individual.n_elem; ++i) {
      // cout << "GetXCoordinate: " << _base_individual(i)(0) << " " << province_id << endl;
        if (_base_individual(i)(0) == province_id) {
            return _base_individual(i)(1);
        }
    }
    return 0.0;
}

double Individual::GetYCoordinate(int province_id) const {
    for (int i = 0; i < _base_individual.n_elem; ++i) {
        if (_base_individual(i)(0) == province_id) {
            return _base_individual(i)(2);
        }
    }
    return 0.0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/