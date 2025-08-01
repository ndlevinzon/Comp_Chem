#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <map>

// Define M_PI if not already defined (for Windows)
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct Atom {
    std::string name;
    std::string resname;
    int resid;
    std::string chain;
    double x, y, z;
    std::string element;
};

struct Residue {
    int resid;
    std::string resname;
    std::vector<Atom> atoms;
};

struct DihedralAngles {
    double phi;
    double psi;
};

class PDBReader {
private:
    std::vector<std::vector<Residue>> frames;
    std::string filename;

    // Helper function to calculate distance between two points
    double distance(const Atom& a, const Atom& b) {
        double dx = a.x - b.x;
        double dy = a.y - b.y;
        double dz = a.z - b.z;
        return sqrt(dx * dx + dy * dy + dz * dz);
    }

    // Helper function to calculate angle between three points
    double angle(const Atom& a, const Atom& b, const Atom& c) {
        double ab = distance(a, b);
        double bc = distance(b, c);
        double ac = distance(a, c);

        if (ab == 0 || bc == 0) return 0.0;

        double cos_angle = (ab * ab + bc * bc - ac * ac) / (2.0 * ab * bc);
        cos_angle = std::max(-1.0, std::min(1.0, cos_angle)); // Clamp to [-1, 1]
        return acos(cos_angle) * 180.0 / M_PI;
    }

    // Helper function to calculate dihedral angle between four points
    double dihedral(const Atom& a, const Atom& b, const Atom& c, const Atom& d) {
        // Calculate vectors
        double v1x = b.x - a.x;
        double v1y = b.y - a.y;
        double v1z = b.z - a.z;

        double v2x = c.x - b.x;
        double v2y = c.y - b.y;
        double v2z = c.z - b.z;

        double v3x = d.x - c.x;
        double v3y = d.y - c.y;
        double v3z = d.z - c.z;

        // Calculate cross products
        double n1x = v1y * v2z - v1z * v2y;
        double n1y = v1z * v2x - v1x * v2z;
        double n1z = v1x * v2y - v1y * v2x;

        double n2x = v2y * v3z - v2z * v3y;
        double n2y = v2z * v3x - v2x * v3z;
        double n2z = v2x * v3y - v2y * v3x;

        // Normalize n1 and n2
        double n1_mag = sqrt(n1x * n1x + n1y * n1y + n1z * n1z);
        double n2_mag = sqrt(n2x * n2x + n2y * n2y + n2z * n2z);

        if (n1_mag == 0 || n2_mag == 0) return 0.0;

        n1x /= n1_mag;
        n1y /= n1_mag;
        n1z /= n1_mag;

        n2x /= n2_mag;
        n2y /= n2_mag;
        n2z /= n2_mag;

        // Calculate dot product
        double dot_product = n1x * n2x + n1y * n2y + n1z * n2z;
        dot_product = std::max(-1.0, std::min(1.0, dot_product)); // Clamp to [-1, 1]

        // Calculate dihedral angle
        double dihedral_angle = acos(dot_product) * 180.0 / M_PI;

        // Determine sign using cross product of n1 and n2
        double sign_x = n1y * n2z - n1z * n2y;
        double sign_y = n1z * n2x - n1x * n2z;
        double sign_z = n1x * n2y - n1y * n2x;

        double sign_dot = sign_x * v2x + sign_y * v2y + sign_z * v2z;

        if (sign_dot < 0) {
            dihedral_angle = -dihedral_angle;
        }

        return dihedral_angle;
    }

    // Helper function to find atom by name in a residue
    Atom* find_atom(Residue& res, const std::string& atom_name) {
        for (auto& atom : res.atoms) {
            if (atom.name == atom_name) {
                return &atom;
            }
        }
        return nullptr;
    }

public:
    PDBReader(const std::string& fname) : filename(fname) {}

    bool readPDB() {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            return false;
        }

        std::string line;
        std::vector<Residue> current_frame;
        std::map<int, Residue> current_residues;

        while (std::getline(file, line)) {
            if (line.substr(0, 6) == "MODEL ") {
                // Start of new frame
                if (!current_residues.empty()) {
                    // Convert map to vector and add to frames
                    for (auto& pair : current_residues) {
                        current_frame.push_back(pair.second);
                    }
                    frames.push_back(current_frame);
                    current_frame.clear();
                    current_residues.clear();
                }
            }
            else if (line.substr(0, 4) == "ATOM" || line.substr(0, 6) == "HETATM") {
                // Parse ATOM/HETATM line
                Atom atom;

                try {
                    atom.name = line.substr(12, 4);
                    atom.resname = line.substr(17, 3);
                    atom.chain = line.substr(21, 1);
                    atom.resid = std::stoi(line.substr(22, 4));
                    atom.x = std::stod(line.substr(30, 8));
                    atom.y = std::stod(line.substr(38, 8));
                    atom.z = std::stod(line.substr(46, 8));
                    atom.element = line.substr(76, 2);

                    // Remove whitespace
                    atom.name.erase(0, atom.name.find_first_not_of(" \t"));
                    atom.name.erase(atom.name.find_last_not_of(" \t") + 1);
                    atom.resname.erase(0, atom.resname.find_first_not_of(" \t"));
                    atom.resname.erase(atom.resname.find_last_not_of(" \t") + 1);
                    atom.chain.erase(0, atom.chain.find_first_not_of(" \t"));
                    atom.chain.erase(atom.chain.find_last_not_of(" \t") + 1);
                    atom.element.erase(0, atom.element.find_first_not_of(" \t"));
                    atom.element.erase(atom.element.find_last_not_of(" \t") + 1);

                    // Add atom to current residue
                    current_residues[atom.resid].resid = atom.resid;
                    current_residues[atom.resid].resname = atom.resname;
                    current_residues[atom.resid].atoms.push_back(atom);
                }
                catch (const std::exception& e) {
                    std::cerr << "Warning: Could not parse line: " << line << std::endl;
                    continue;
                }
            }
            else if (line.substr(0, 6) == "ENDMDL") {
                // End of current frame
                if (!current_residues.empty()) {
                    // Convert map to vector and add to frames
                    for (auto& pair : current_residues) {
                        current_frame.push_back(pair.second);
                    }
                    frames.push_back(current_frame);
                    current_frame.clear();
                    current_residues.clear();
                }
            }
        }

        // Handle last frame if no ENDMDL
        if (!current_residues.empty()) {
            for (auto& pair : current_residues) {
                current_frame.push_back(pair.second);
            }
            frames.push_back(current_frame);
        }

        file.close();
        std::cout << "Read " << frames.size() << " frames from PDB file." << std::endl;
        return true;
    }

    DihedralAngles calculateDihedrals(int target_resid, int frame_idx) {
        DihedralAngles angles = { 0.0, 0.0 };

        if (frame_idx >= static_cast<int>(frames.size())) {
            std::cerr << "Error: Frame index out of range." << std::endl;
            return angles;
        }

        // Find target residue in current frame
        Residue* target_res = nullptr;
        for (auto& res : frames[frame_idx]) {
            if (res.resid == target_resid) {
                target_res = &res;
                break;
            }
        }

        if (!target_res) {
            std::cerr << "Error: Residue " << target_resid << " not found in frame " << frame_idx << std::endl;
            return angles;
        }

        // Find atoms needed for phi calculation: C(i-1) - N(i) - CA(i) - C(i)
        Atom* C_prev = nullptr, * N = nullptr, * CA = nullptr, * C = nullptr;

        // Find N, CA, C in current residue
        N = find_atom(*target_res, "N");
        CA = find_atom(*target_res, "CA");
        C = find_atom(*target_res, "C");

        if (!N || !CA || !C) {
            std::cerr << "Error: Missing backbone atoms in residue " << target_resid << std::endl;
            return angles;
        }

        // Find C from previous residue
        for (auto& res : frames[frame_idx]) {
            if (res.resid == target_resid - 1) {
                C_prev = find_atom(res, "C");
                break;
            }
        }

        // Calculate phi if we have all required atoms
        if (C_prev && N && CA && C) {
            angles.phi = dihedral(*C_prev, *N, *CA, *C);
        }

        // Find atoms needed for psi calculation: N(i) - CA(i) - C(i) - N(i+1)
        Atom* N_next = nullptr;

        // Find N from next residue
        for (auto& res : frames[frame_idx]) {
            if (res.resid == target_resid + 1) {
                N_next = find_atom(res, "N");
                break;
            }
        }

        // Calculate psi if we have all required atoms
        if (N && CA && C && N_next) {
            angles.psi = dihedral(*N, *CA, *C, *N_next);
        }

        return angles;
    }

    int getNumFrames() const {
        return static_cast<int>(frames.size());
    }

    void writeCSV(const std::string& output_filename, int target_resid, double timestep_ps = 100.0) {
        std::ofstream csv_file(output_filename);
        if (!csv_file.is_open()) {
            std::cerr << "Error: Could not create output file " << output_filename << std::endl;
            return;
        }

        // Write header
        csv_file << "frame,time_ps,phi_deg,psi_deg\n";

        // Write data for each frame
        for (int frame = 0; frame < static_cast<int>(frames.size()); ++frame) {
            DihedralAngles angles = calculateDihedrals(target_resid, frame);
            double time_ps = frame * timestep_ps;

            csv_file << frame << ","
                << std::fixed << std::setprecision(2) << time_ps << ","
                << std::fixed << std::setprecision(3) << angles.phi << ","
                << std::fixed << std::setprecision(3) << angles.psi << "\n";

            if (frame % 100 == 0) {
                std::cout << "Processed frame " << frame << "/" << frames.size() - 1 << std::endl;
            }
        }

        csv_file.close();
        std::cout << "Results written to " << output_filename << std::endl;
    }
};

int main() {
    std::string pdb_filename;
    int target_residue;
    double timestep_ps;

    // Get input from user
    std::cout << "Enter PDB filename: ";
    std::getline(std::cin, pdb_filename);

    std::cout << "Enter target residue number: ";
    std::cin >> target_residue;

    std::cout << "Enter timestep in picoseconds (default 100.0): ";
    std::cin >> timestep_ps;

    // Create PDB reader and process file
    PDBReader reader(pdb_filename);

    if (!reader.readPDB()) {
        std::cerr << "Failed to read PDB file." << std::endl;
        return 1;
    }

    // Generate output filename
    std::string output_filename = pdb_filename;
    size_t dot_pos = output_filename.find_last_of('.');
    if (dot_pos != std::string::npos) {
        output_filename = output_filename.substr(0, dot_pos);
    }
    output_filename += "_res" + std::to_string(target_residue) + "_phi_psi.csv";

    // Calculate dihedrals and write to CSV
    reader.writeCSV(output_filename, target_residue, timestep_ps);

    std::cout << "Analysis complete!" << std::endl;
    return 0;
}
