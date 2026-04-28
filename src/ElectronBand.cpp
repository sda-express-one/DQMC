#include "../include/ElectronBand.hpp"

ParabolicBandSimple::ParabolicBandSimple(double effective_mass) : _effective_mass(effective_mass) {};

ParabolicBandAnisotropic::ParabolicBandAnisotropic(double mx_el, double my_el, double mz_el) : _mx_el(mx_el), _my_el(my_el), _mz_el(mz_el) {};

ParabolicBandLK::ParabolicBandLK(double A_LK, double B_LK, double C_LK) : _A_LK(A_LK), _B_LK(B_LK), _C_LK(C_LK) {};

FullBand::FullBand(const std::string& filename) {
    // Load the full band structure from the file
    double kx = 0, ky = 0, kz = 0, energy = 0;
    

    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file: " << filename << std::endl;
        //throw std::runtime_error("Could not open file: " + filename);
    }

    std::string line;

    std::string cell_type = "None";
    int num_kpoints = 1;
    int num_bands = 1;
    int lineNum = 0;
    bool acquiring_data = false;
    int i = 0;
    
    while(std::getline(file, line)){
        if(line.empty()){continue;}

        // remove everything after # if present
        size_t commentPosition = line.find("#");
        if(commentPosition != std::string::npos){
            line = line.substr(0, commentPosition);
        }
        
        if(acquiring_data){
            // retrieve band values
            std::stringstream ss(line);
            ss >> _kpoints_table[3*i] >> _kpoints_table[3*i+1] >> _kpoints_table[3*i+2];
            for(int j = 0; j < _num_bands; ++j){
                ss >> _energy_table[_num_bands*i+j];
            }
            ++i;
        }
        else{
            // first uncommented line: get cell type
            if(lineNum == 0){
                std::stringstream ss(line);
                ss >> cell_type;
                cell_type = _cell_type;
            }
            // second uncommented line: get 2 values
            if(lineNum == 1){
                std::stringstream ss(line);
                ss >> num_bands >> num_kpoints;
                _num_bands = num_bands;
                _num_kpoints = num_kpoints;
                _bands_type = new int[num_bands];
                _kpoints_table = new double[num_kpoints * 3];
                _energy_table = new double[num_bands * num_kpoints];
                ++lineNum;
            }
            // third to fifth uncommented line: get lattice vectors
            else if(lineNum > 1 && lineNum < 5){
                std::stringstream ss(line);
                ss >> _unit_cell(lineNum-2,0) >> _unit_cell(lineNum-2,1) >> _unit_cell(lineNum-2,2);
                ++lineNum;
            }
            // sixth uncommented line: get bands type (valence or conduction)
            else if(lineNum == 5){
                std::stringstream ss(line);
                for(int i = 0; i < num_bands; ++i){
                    ss >> _bands_type[i];
                }
                ++lineNum;
            }
            else{
                std::stringstream ss(line);
                std::string header_check;
                ss >> header_check;
                if(header_check == "HEADER END"){
                    acquiring_data = true;
                };
            }
        }
    }
    file.close();

    double energy_extremum = 0;
    int idx = 0;

    for(int band=0; band<_num_bands; ++band){
        
        if(_bands_type[0] == 1){
            // If the bands are conduction bands, shift energies so that the bottom of the conduction band is at 0
            for (int i = 0; i < _num_kpoints * _num_bands; ++i) {
                idx = i * 4 + 3;
                _energy_table[idx] -= energy_extremum;
            }
        }
        else if(_bands_type[0] == 0){
            // If the bands are valence bands, invert the energies so that the top of the valence band is at 0
            for (int i = 0; i < _num_kpoints * _num_bands; ++i) {
                idx = i * 4 + 3;
                _energy_table[idx] = energy_extremum - _energy_table[idx];
            }
        }
    }
};

void FullBand::buildUnitCellRec() {
    Eigen::Vector3d basis_vector_rep;
    double volume = 1;
    Eigen::Vector3d a1 = _unit_cell.col(0), a2 = _unit_cell.col(1), a3 = _unit_cell.col(2);
    volume = a1.dot(a2.cross(a3));
    _unit_cell_rec.col(0) = 2*M_PI/volume*a2.cross(a3);
    _unit_cell_rec.col(1) = 2*M_PI/volume*a3.cross(a1);
    _unit_cell_rec.col(2) = 2*M_PI/volume*a1.cross(a2);
};

int FullBand::findClosestKPoint(double kx, double ky, double kz) const {
    // placeholder
    return 0;
};