//-------------------------------------------------------
//
//  Unit Cell creator
//
//  Takes a primitive unit cell and replicates it, 
//  creating the neighbourlist and populating atomic 
//  properties for input into vampire
//
//  (C) R.F.L.Evans 22/04/2015
//
//
//-------------------------------------------------------
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

class uc_atom_t{

public:

  double cx;
  double cy;
  double cz;

  int material;
  int hc;
  int lc;

};

class material_t{

public:

  double mu_s;
  double alpha;
  double Ku;
  double Sx;
  double Sy;
  double Sz;
  std::string name;
  std::string element;


};

class nn_t{

public:
  int i;
  int j;
  int dx;
  int dy;
  int dz;
  double Jij = 0.0;
  double Dx = 0.0;
  double Dy = 0.0;
  double Dz = 0.0;
  double rij;
};

int main(){

  // system constants
  //unit cell sizes
  double unit_cell_size[3]={7.972, 6.904, 16.229};

  const int num_materials=6;
  // exchange constants
  std::vector<std::vector<double> > exchange_constants;
  exchange_constants.resize(num_materials);
  for(int m=0;m<num_materials;m++) exchange_constants.at(m).resize(num_materials, 0.0);
  std::vector<std::vector<double> > exchange_constants_nnn;
  exchange_constants_nnn.resize(num_materials);
  for(int m=0;m<num_materials;m++) exchange_constants_nnn.at(m).resize(num_materials, 0.0);

  // material parameters
  std::vector<material_t> materials(2);

    
  materials.at(0).mu_s=2.05; // mu_B's
  materials.at(0).alpha=0.1;
  materials.at(0).name="Fe1,2";
  materials.at(0).element="Fe";
  materials.at(0).Sx=0.0;
  materials.at(0).Sy=0.0;
  materials.at(0).Sz=1.0;

  materials.at(1).mu_s=1.35;
  materials.at(1).alpha=0.1;
  materials.at(1).name="Fe3";
  materials.at(1).element="Fe";
  materials.at(1).Sx=0.0;
  materials.at(1).Sy=0.0;
  materials.at(1).Sz=1.0;
  // create atoms in unit cell
  std::vector<uc_atom_t> unit_cell(24);

   double x1 = 0.125;
   double x2 = 0.375;
   double x3 = 0.625;
   double x4 = 0.875;

   double y1 = 0.0833574739281576;
   double y2 = 0.250072421784473;
   double y3 = 0.416642526071842;
   double y4 = 0.583357473928158;
   double y5 = 0.750072421784473;
   double y6 = 0.916642526071842;

   double z1 = 0.176967157557459; //type 0
   double z2 = 0.249984595477232 + 0.00544375; //type 2 //noncentroshift 0.0871 A / 16.0 A -> 0.00544375
   double z3 = 0.323063651488077; //type 1
   double z4 = 0.676936348511923; //type 3
   double z5 = 0.750015404522768 + 0.00544375; //type 5 //noncentroshift 0.0871 A / 16.0 A -> 0.00544375
   double z6 = 0.823032842442541; //type 4

  unit_cell.at(0).cx=x1;
  unit_cell.at(0).cy=y2;
  unit_cell.at(0).cz=z5;
  unit_cell.at(0).material=1;
  unit_cell.at(0).hc=4;
  unit_cell.at(0).lc=6;
  
  unit_cell.at(1).cx=x2;
  unit_cell.at(1).cy=y1;
  unit_cell.at(1).cz=z1;
  unit_cell.at(1).material=0;
  unit_cell.at(1).hc=1;
  unit_cell.at(1).lc=1;
  
  unit_cell.at(2).cx=x2;
  unit_cell.at(2).cy=y3;
  unit_cell.at(2).cz=z2;
  unit_cell.at(2).material=1;
  unit_cell.at(2).hc=1;
  unit_cell.at(2).lc=3;
  
  unit_cell.at(3).cx=x2;
  unit_cell.at(3).cy=y1;
  unit_cell.at(3).cz=z3;
  unit_cell.at(3).material=0;
  unit_cell.at(3).hc=1;
  unit_cell.at(3).lc=2;
  
  unit_cell.at(4).cx=x2;
  unit_cell.at(4).cy=y1;
  unit_cell.at(4).cz=z4;
  unit_cell.at(4).material=0;
  unit_cell.at(4).hc=1;
  unit_cell.at(4).lc=4;

  unit_cell.at(5).cx=x3;
  unit_cell.at(5).cy=y2;
  unit_cell.at(5).cz=z5;
  unit_cell.at(5).material=1;
  unit_cell.at(5).hc=4;
  unit_cell.at(5).lc=6;

  unit_cell.at(6).cx=x2;
  unit_cell.at(6).cy=y1;
  unit_cell.at(6).cz=z6;
  unit_cell.at(6).material=0;
  unit_cell.at(6).hc=1;
  unit_cell.at(6).lc=5;

  unit_cell.at(7).cx=x1;
  unit_cell.at(7).cy=y4;
  unit_cell.at(7).cz=z1;
  unit_cell.at(7).material=0;
  unit_cell.at(7).hc=1;
  unit_cell.at(7).lc=1;

  unit_cell.at(8).cx=x1;
  unit_cell.at(8).cy=y6;
  unit_cell.at(8).cz=z2;
  unit_cell.at(8).material=1;
  unit_cell.at(8).hc=1;
  unit_cell.at(8).lc=3;

  unit_cell.at(9).cx=x1;
  unit_cell.at(9).cy=y4;
  unit_cell.at(9).cz=z3;
  unit_cell.at(9).material=0;
  unit_cell.at(9).hc=1;
  unit_cell.at(9).lc=2;

  unit_cell.at(10).cx=x1;
  unit_cell.at(10).cy=y4;
  unit_cell.at(10).cz=z4;
  unit_cell.at(10).material=0;
  unit_cell.at(10).hc=1;
  unit_cell.at(10).lc=4;

  unit_cell.at(11).cx=x2;
  unit_cell.at(11).cy=y5;
  unit_cell.at(11).cz=z5;
  unit_cell.at(11).material=1;
  unit_cell.at(11).hc=1;
  unit_cell.at(11).lc=6;

  unit_cell.at(12).cx=x1;
  unit_cell.at(12).cy=y4;
  unit_cell.at(12).cz=z6;
  unit_cell.at(12).material=0;
  unit_cell.at(12).hc=1;
  unit_cell.at(12).lc=5;

  unit_cell.at(13).cx=x4;
  unit_cell.at(13).cy=y1;
  unit_cell.at(13).cz=z1;
  unit_cell.at(13).material=0;
  unit_cell.at(13).hc=1;
  unit_cell.at(13).lc=1;

  unit_cell.at(14).cx=x4;
  unit_cell.at(14).cy=y3;
  unit_cell.at(14).cz=z2;
  unit_cell.at(14).material=1;
  unit_cell.at(14).hc=1;
  unit_cell.at(14).lc=3;

  unit_cell.at(15).cx=x4;
  unit_cell.at(15).cy=y1;
  unit_cell.at(15).cz=z3;
  unit_cell.at(15).material=0;
  unit_cell.at(15).hc=1;
  unit_cell.at(15).lc=2;

  unit_cell.at(16).cx=x4;
  unit_cell.at(16).cy=y1;
  unit_cell.at(16).cz=z4;
  unit_cell.at(16).material=0;
  unit_cell.at(16).hc=1;
  unit_cell.at(16).lc=4;

  unit_cell.at(17).cx=x4;
  unit_cell.at(17).cy=y1;
  unit_cell.at(17).cz=z6;
  unit_cell.at(17).material=0;
  unit_cell.at(17).hc=1;
  unit_cell.at(17).lc=5;

  unit_cell.at(18).cx=x3;
  unit_cell.at(18).cy=y4;
  unit_cell.at(18).cz=z1;
  unit_cell.at(18).material=0;
  unit_cell.at(18).hc=1;
  unit_cell.at(18).lc=1;

  unit_cell.at(19).cx=x3;
  unit_cell.at(19).cy=y6;
  unit_cell.at(19).cz=z2;
  unit_cell.at(19).material=1;
  unit_cell.at(19).hc=1;
  unit_cell.at(19).lc=3;

  unit_cell.at(20).cx=x3;
  unit_cell.at(20).cy=y4;
  unit_cell.at(20).cz=z3;
  unit_cell.at(20).material=0;
  unit_cell.at(20).hc=1;
  unit_cell.at(20).lc=2;

  unit_cell.at(21).cx=x3;
  unit_cell.at(21).cy=y4;
  unit_cell.at(21).cz=z4;
  unit_cell.at(21).material=0;
  unit_cell.at(21).hc=1;
  unit_cell.at(21).lc=4;
  
  unit_cell.at(22).cx=x4;
  unit_cell.at(22).cy=y5;
  unit_cell.at(22).cz=z5;
  unit_cell.at(22).material=1;
  unit_cell.at(22).hc=1;
  unit_cell.at(22).lc=6;

  unit_cell.at(23).cx=x3;
  unit_cell.at(23).cy=y4;
  unit_cell.at(23).cz=z6;
  unit_cell.at(23).material=0;
  unit_cell.at(23).hc=1;
  unit_cell.at(23).lc=5;

   //unit_cell.cutoff_radius = 0.15;
  double cell_x_norm = 1.0 ;//0.49121942;
  double cell_y_norm = 1.0 ;//0.4254113;
  double cell_z_norm = 1.0;
   //normalise distance to tetragonal cell

    //GGA constants, DFT (e.g., needs 2x factor)
  // double J12 = 35.52;
  // double J13 = 10.42;
  // double J11 = 3.25;
  // double J33 = 0.64;
  // double J12_nnn = 2.18;
  // double J13_nnn = 0.77;
  // double J24_nnn = 0.29;

    //LDA
  // double J12 = 30.72;
  // double J13 = 12.19;
  // double J11 = 4.68;
  // double J33 = -0.04;
  // double J12_nnn = 1.53;
  // double J13_nnn = 1.30;
  // double J24_nnn = 0.68;

  //   //Li et al., (LDA)
  // double J12 = 57.18; 
  // double J13 = 17.02;
  // double J11 = -0.92;
  // double J33 = 2.66; //the Te mediated DMI one. DMI: 1 meV (10.1038/s41467-024-48799-9)
  // double J12_nnn = -1.55;
  // double J13_nnn = -1.29;
  // double J24_nnn = 3.43;

  double DMI_exchange_ratio = 0.00;
    //Shreyas 
  double J12 = 85.3134; //<-0.000,  0.000,  2.405> ( 0.0046  0.0011 -0.0000)
  double J13 = 14.7714; // ( 2.044,  1.177, -1.201) (-0.1940  0.3483 -0.0805)
  double J11 = -0.92;
  double J33 = 2.66;
  double J12_nnn = -1.55;
  double J13_nnn = -1.29;
  double J24_nnn = 3.43;

  exchange_constants[0][0] = J11;
  exchange_constants[0][1] = J12;
  exchange_constants[0][2] = J13;

  exchange_constants[1][0] = J12;
  exchange_constants[1][1] = J11;
  exchange_constants[1][2] = J13;

  exchange_constants[2][0] = J13;
  exchange_constants[2][1] = J13;
  exchange_constants[2][2] = J33;

  exchange_constants[3][3] = J11;
  exchange_constants[3][4] = J12;
  exchange_constants[3][5] = J13;

  exchange_constants[4][3] = J12;
  exchange_constants[4][4] = J11;
  exchange_constants[4][5] = J13;

  exchange_constants[5][3] = J13;
  exchange_constants[5][4] = J13;
  exchange_constants[5][5] = J33;
  
  exchange_constants_nnn[0][1] = J12_nnn;
  exchange_constants_nnn[0][2] = J13_nnn;
  exchange_constants_nnn[0][4] = J24_nnn;

  exchange_constants_nnn[1][0] = J12_nnn;
  exchange_constants_nnn[1][2] = J13_nnn;
  exchange_constants_nnn[1][3] = J24_nnn;

  exchange_constants_nnn[2][0] = J13_nnn;
  exchange_constants_nnn[2][1] = J13_nnn;

  exchange_constants_nnn[3][4] = J12_nnn;
  exchange_constants_nnn[3][5] = J13_nnn;
  exchange_constants_nnn[3][1] = J24_nnn;

  exchange_constants_nnn[4][3] = J12_nnn;
  exchange_constants_nnn[4][5] = J13_nnn;
  exchange_constants_nnn[4][0] = J24_nnn;

  exchange_constants_nnn[5][3] = J13_nnn;
  exchange_constants_nnn[5][4] = J13_nnn;


  for(int i = 0; i < unit_cell.size(); i++) {
    unit_cell.at(i).cx *= cell_x_norm*unit_cell_size[0];
    unit_cell.at(i).cy *= cell_y_norm*unit_cell_size[1];
    unit_cell.at(i).cz *= cell_z_norm*unit_cell_size[2];
    unit_cell.at(i).lc--; //reduce by 1 for array math
  }
  // store vector of unit cells
  std::vector< std::vector < std::vector < std::vector<uc_atom_t > > > >crystal;
  
  crystal.resize(3);
  for(int i=0;i<3;i++){
    crystal.at(i).resize(3);
    for(int j=0;j<3;j++){
      crystal.at(i).at(j).resize(3);
      for(int k=0;k<3;k++){
	    crystal.at(i).at(j).at(k).resize(unit_cell.size());
      }
    }
  }
    
  // replicate unit cell
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
	for(int a=0;a<unit_cell.size();a++){
	  crystal.at(i).at(j).at(k).at(a).cx=unit_cell.at(a).cx + unit_cell_size[0]*cell_x_norm*double(i);
    crystal.at(i).at(j).at(k).at(a).cy=unit_cell.at(a).cy + unit_cell_size[1]*cell_y_norm*double(j);
	  crystal.at(i).at(j).at(k).at(a).cz=unit_cell.at(a).cz + unit_cell_size[2]*cell_z_norm*double(k);
	  
    crystal.at(i).at(j).at(k).at(a).material=unit_cell.at(a).material;
	  crystal.at(i).at(j).at(k).at(a).hc=unit_cell.at(a).hc+2*k-2;
    crystal.at(i).at(j).at(k).at(a).lc=unit_cell.at(a).lc;	  
	  }
      }
    }
  }

  // create neighbour list
  double nn_range= 5.75;
  nn_range *= nn_range;
  std::vector<nn_t> nn_list;
  double dmi_ref_vector[3] = {0.0,0.0,1.0};
  // loop over all atoms in unit cell
  for(int ai=0;ai<unit_cell.size();ai++){
    double icx=crystal[1][1][1].at(ai).cx;
    double icy=crystal[1][1][1].at(ai).cy;
    double icz=crystal[1][1][1].at(ai).cz;
    int imat=crystal[1][1][1].at(ai).material;
    int ihc=crystal[1][1][1].at(ai).hc;
    int ilc=crystal[1][1][1].at(ai).lc;

    // loop over all other atoms
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
	for(int k=0;k<3;k++){
	  for(int aj=0;aj<unit_cell.size();aj++){
	    double jcx=crystal[i][j][k].at(aj).cx;
	    double jcy=crystal[i][j][k].at(aj).cy;
	    double jcz=crystal[i][j][k].at(aj).cz;
	    int jmat=crystal[i][j][k].at(aj).material;
	    int jhc=crystal[i][j][k].at(aj).hc;
	    int jlc=crystal[i][j][k].at(aj).lc;
	    double range_sq=(jcx-icx)*(jcx-icx)+(jcy-icy)*(jcy-icy)+(jcz-icz)*(jcz-icz);
	    bool same_atom=(ai==aj && i==1 && j==1 && k==1);
	    if(range_sq<=nn_range && same_atom==false){
	      nn_t temp;
	      temp.i=ai;
	      temp.j=aj;
	      temp.dx=i-1;
	      temp.dy=j-1;
	      temp.dz=k-1;

        double r_inv = 1.0/sqrt(range_sq);
        
        double unit_ij[3] = {(icx-jcx)*r_inv, (icy-jcy)*r_inv, (icz-jcz)*r_inv };
        double uij_cross_z[3] = {unit_ij[2]*dmi_ref_vector[1]-unit_ij[1]*dmi_ref_vector[2], \
                                 unit_ij[0]*dmi_ref_vector[2]-unit_ij[2]*dmi_ref_vector[0], \
                                 unit_ij[1]*dmi_ref_vector[0]-unit_ij[0]*dmi_ref_vector[1]};
	      if(range_sq <= (3.99*3.99)) {
          temp.Jij = exchange_constants.at(ilc).at(jlc);
          if(imat != jmat) {
            temp.Dx = DMI_exchange_ratio* exchange_constants.at(ilc).at(jlc)*uij_cross_z[0];
            temp.Dy = DMI_exchange_ratio* exchange_constants.at(ilc).at(jlc)*uij_cross_z[1];
            temp.Dz = DMI_exchange_ratio* exchange_constants.at(ilc).at(jlc)*uij_cross_z[2];
          } else if (2 == ilc == jlc || 5 == ilc == jlc) {
            // temp.Dx = DMI_exchange_ratio* exchange_constants.at(ilc).at(jlc)*uij_cross_z[0];
            // temp.Dy = DMI_exchange_ratio* exchange_constants.at(ilc).at(jlc)*uij_cross_z[1];
            temp.Dz = DMI_exchange_ratio* exchange_constants.at(ilc).at(jlc);
          }
        }
        else {
          temp.Jij=exchange_constants_nnn.at(ilc).at(jlc);
          // temp.Dx = DMI_exchange_ratio* exchange_constants_nnn.at(ilc).at(jlc)*uij_cross_z[0];
          // temp.Dy = DMI_exchange_ratio* exchange_constants_nnn.at(ilc).at(jlc)*uij_cross_z[1];
          // temp.Dz = DMI_exchange_ratio* exchange_constants_nnn.at(ilc).at(jlc)*uij_cross_z[2];
        }
        if(temp.Jij == 0.0) continue;
        temp.Jij *= 2.0;///materials[imat].mu_s/materials[jmat].mu_s;
        temp.Dx  *= 2.0;
        temp.Dy  *= 2.0;
        temp.Dz  *= 2.0;

        temp.rij = sqrt(range_sq);
	      nn_list.push_back(temp);
	      //std::cout << ai << "\t" << aj << "\t" << i-1 << "\t" << j-1 << "\t" << k-1 << "\t" << sqrt(range_sq) << std::endl;
	    } 
	  }
	}
      }
    }

  }


  // output to files
  // declare outfile file stream
  std::ofstream ucf_file;
  // open it (file_name)
  ucf_file.open ("fe3gate2-lda-Li-tetragonal-ncsDMI-10percent.ucf");

  ucf_file << "# Unit cell size:" << std::endl;
  ucf_file << unit_cell_size[0]*cell_x_norm << "\t" << unit_cell_size[1]*cell_y_norm << "\t" << unit_cell_size[2] << std::endl;
  ucf_file << "# Unit cell vectors: " << std::endl;
  ucf_file << "1.0 0.0 0.0 " << std::endl;
  ucf_file << "0.0 1.0 0.0 " << std::endl;
  ucf_file << "0.0 0.0 1.0 " << std::endl;
  ucf_file << "# Atoms num, id cx cy cz mat lc hc " << std::endl;
  ucf_file << unit_cell.size() << std::endl;
  // loop over all atoms
  for(int atom=0; atom<unit_cell.size(); atom++){
    ucf_file << atom << "\t";
    ucf_file << unit_cell.at(atom).cx/unit_cell_size[0]/cell_x_norm << "\t";
    ucf_file << unit_cell.at(atom).cy/unit_cell_size[1]/cell_y_norm << "\t";
    ucf_file << unit_cell.at(atom).cz/unit_cell_size[2]/cell_z_norm << "\t";
    ucf_file << unit_cell.at(atom).material << "\t";
    ucf_file << unit_cell.at(atom).lc << "\t";
    ucf_file << unit_cell.at(atom).hc << std::endl;
  }
  ucf_file << "#Interactions n exctype, id i j dx dy   dz        Jij"<< std::endl;
  ucf_file << nn_list.size() << "\t" << "normalised-isotropic" << std::endl;
  // loop over all interactions

  for(unsigned int nn=0; nn < nn_list.size(); nn++){

    ucf_file << nn << "\t";
    ucf_file << nn_list[nn].i << "\t";
    ucf_file << nn_list[nn].j << "\t"; 
    ucf_file << nn_list[nn].dx << "\t";
    ucf_file << nn_list[nn].dy << "\t"; 
    ucf_file << nn_list[nn].dz << "\t"; 
    ucf_file <<  nn_list[nn].Jij << "\t" << nn_list[nn].Dz  << "\t" << -nn_list[nn].Dy << "\t" \
             << -nn_list[nn].Dz  << "\t" << nn_list[nn].Jij << "\t" <<  nn_list[nn].Dx << "\t" \
             << -nn_list[nn].Dy  << "\t" <<-nn_list[nn].Dx  << "\t" <<  nn_list[nn].Jij << "\t" << \
    " #" << nn_list[nn].rij << std::endl;

  }
  
  // material file
  std::ofstream mat_file;
  // open it (file_name)
  mat_file.open ("MAT.material");
  mat_file << "#================================================" << std::endl;
  mat_file << "# Generated material file for input into vampire" << std::endl;
  mat_file << "#================================================" << std::endl;
  mat_file << "#" << std::endl;
  mat_file << "# File timestamp: " << std::endl;
  mat_file << "#" << std::endl;
  mat_file << "#------------------------------------------------" << std::endl;
  mat_file << "# Number of Materials" << std::endl;
  mat_file << "#------------------------------------------------" << std::endl;
  mat_file << "material:num-materials=" << num_materials << std::endl;
  mat_file << "#------------------------------------------------" << std::endl;
  
  // Loop over all materials
  for(int m=0;m<materials.size();m++){
    mat_file << "# Material " << m+1 << " (" << materials.at(m).name << ")" << std::endl;
    mat_file << "#------------------------------------------------" << std::endl;
    mat_file << "material[" << m+1 << "]:material-name=" << materials.at(m).name << std::endl;
    mat_file << "material[" << m+1 << "]:damping-constant=" << materials.at(m).alpha << std::endl;
    mat_file << "material[" << m+1 << "]:atomic-spin-moment="<< materials.at(m).mu_s << " !muB" << std::endl;
    mat_file << "material[" << m+1 << "]:uniaxial-anisotropy-constant=" << materials.at(m).Ku << std::endl;
    //mat_file << "material[" << m << "]:uniaxial-anisotropy-direction=" << 0 << "," << 1 << ","<< 0 << std::endl;
    mat_file << "material[" << m+1 << "]:material-element="<< materials.at(m).element << ""<< std::endl;
    mat_file << "material[" << m+1 << "]:initial-spin-direction="<< materials.at(m).Sx << "," << materials.at(m).Sy << "," << materials.at(m).Sz << std::endl;
    mat_file << "#------------------------------------------------" << std::endl;
  }
  return 0;
}
