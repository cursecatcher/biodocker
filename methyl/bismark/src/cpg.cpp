#include <algorithm>
#include <iostream>
#include <fstream>
#include <map>
#include <list>
#include <chrono>
#include <vector>

struct Coordinate {
    const std::string chromosome; 
    const unsigned int position; 

    Coordinate(std::string& chr, unsigned int pos) 
        : chromosome(chr), position(pos) {
    }

    inline std::string to_string() {
        return chromosome + ":" + std::to_string(position);
    }
};

class CpGRegion {
    unsigned int count;
public:
    const std::string chromosome; 
    const unsigned int start, end; 

    CpGRegion(std::string& chr, unsigned int start, unsigned int end) 
        : chromosome(chr), start(start), end(end), count(0) {
    }

    inline void increment() { 
        ++count; 
    }

    inline int get_count() { 
        return count; 
    }
};



class CpGRegions {
    std::map<std::string, CpGRegion&> regions; 
    std::list<CpGRegion> region_list; 
public: 
    CpGRegions() {
    }

    void add(std::string& chr, unsigned start, unsigned end) {
        region_list.emplace_back(chr, start, end); 
        CpGRegion& region = region_list.back(); 

        for (unsigned pos = start; pos <= end; ++pos) {
            regions.emplace(
                std::piecewise_construct, 
                std::forward_as_tuple(Coordinate(chr, pos).to_string()), 
                std::forward_as_tuple(region)
            );
        }
    }

    inline CpGRegion* match(Coordinate coord) {
        auto it = regions.find(coord.to_string()); 
        return it == regions.end() ? nullptr : &it->second; 
    }

    inline std::list<CpGRegion>::iterator begin() {
        return region_list.begin(); 
    }

    inline std::list<CpGRegion>::iterator end() {
        return region_list.end(); 
    }
}; 



int main(int argc, char** argv) {
    std::string cpg_region_file(argv[1]);
    std::string methylation_file(argv[2]); 
    std::string cpg_count_per_region = "cpg_count_per_region.bed"; 
    std::string off_targets = "off_targets.bed";

    
    std::string tmp; 
    std::string chr, cpg_id; 
    unsigned start, end; 

    float tot = 0, n_regions = 0; 

    CpGRegions my_regions;

    std::cout << "Creating hash table for CpG regions..." << std::endl; 
    std::ifstream cpg_file(cpg_region_file); 
    //std::getline(cpg_file, tmp);
    auto start_t = std::chrono::steady_clock::now();

    for (; cpg_file >> chr >> start >> end >> cpg_id; ++n_regions, tot += end - start) 
        my_regions.add(chr, start, end);

    cpg_file.close(); 

    auto end_t = std::chrono::steady_clock::now(); 
    std::cout 
        << " Completed in " 
        << std::chrono::duration_cast<std::chrono::seconds>(end_t-start_t).count() << "seconds\n"
        << "#regions: " << n_regions << "\nAvg size: " << tot / n_regions << std::endl;
    

    std::ifstream meth_file(methylation_file); 
    std::getline(meth_file, tmp); 
    std::string read_id; 
    char strand, status; 

    std::ofstream off_targets_file(off_targets); 

    std::cout << "Matching methylation coordinates..." << std::endl; 
    start_t = std::chrono::steady_clock::now(); 
    for (; meth_file >>  read_id >> strand >> chr >> start >> status;) {
        CpGRegion* region = my_regions.match(Coordinate(chr, start)); 
        if (region != nullptr)
            region->increment();
        else 
            off_targets_file << chr << "\t" << start << "\t" << strand << "\t" << status << "\n"; 
    }
    end_t = std::chrono::steady_clock::now(); 
    std::cout 
        << "Completed in " 
        << std::chrono::duration_cast<std::chrono::seconds>(end_t-start_t).count() << "seconds" << std::endl; 

    std::ofstream cpg_count_per_region_file(cpg_count_per_region); 
    std::vector<CpGRegion> v_region{std::begin(my_regions), std::end(my_regions)}; 
    /*
    std::sort(v_region.begin(), v_region.end(), [](const CpGRegion& x, const CpGRegion& y) -> bool {
        return x.chromosome <= y.chromosome; 
    }); */

    for (auto it = v_region.begin(); it != v_region.end(); ++it) {
        cpg_count_per_region_file 
            << it->chromosome  << "\t" 
            << it->start << "\t" 
            << it->end << "\t" 
            << it->get_count() << "\n"; 
    }

    unsigned count;
    for (auto it = my_regions.begin(); it != my_regions.end(); ++it) { 
        if ((count = it->get_count()) > 0) {
            std::cout << it->chromosome << ":" << it->start << "_" << it->end << " -> " << count << std::endl; 
        }
    }
}