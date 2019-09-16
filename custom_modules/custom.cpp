/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include <sstream>
#include <string>

#include "./custom.h"

// declare cell definitions here 

void create_cell_types(void)
{
	static int idxDeath = microenvironment.find_density_index("death_signal");
	SeedRandom(parameters.ints("random_seed")); // or specify a seed here 

	initialize_default_cell_definition();

	// Name the default cell type 
	cell_defaults.type = 0;
	cell_defaults.name = "my cell";
	cell_defaults.phenotype.mechanics.cell_cell_adhesion_strength = 0.0;

	// make sure the defaults are self-consistent.  (rwh: do this BEFORE others)
	cell_defaults.phenotype.secretion.sync_to_microenvironment(&microenvironment);
	cell_defaults.phenotype.molecular.sync_to_microenvironment(&microenvironment);

	// set default cell cycle model 
	// cell_defaults.functions.cycle_model = flow_cytometry_separated_cycle_model; 

	// set default_cell_functions; 
	cell_defaults.functions.update_phenotype = death_function;

	// release particles at death 
	cell_defaults.phenotype.molecular.fraction_released_at_death[idxDeath] = 1.0;

	cell_defaults.phenotype.secretion.uptake_rates[idxDeath] = parameters.doubles("signal_internalization_rate");

	// needed for a 2-D simulation: 

	cell_defaults.functions.set_orientation = up_orientation;
	cell_defaults.phenotype.geometry.polarity = 1.0;
	cell_defaults.phenotype.motility.restrict_to_2D = true;


	cell_defaults.phenotype.sync_to_functions(cell_defaults.functions);

	// set the rate terms in the default phenotype 

	// no death (we'll lyse cell ourself in the custom 'death_function')
	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index("Apoptosis");
	int necrosis_model_index = cell_defaults.phenotype.death.find_death_model_index("Necrosis");
	cell_defaults.phenotype.death.rates[apoptosis_model_index] = 0.0;
	cell_defaults.phenotype.death.rates[necrosis_model_index] = 0.0;

	// no birth (turn off proliferation)
	int cycle_start_index = live.find_phase_index(PhysiCell_constants::live);
	int cycle_end_index = live.find_phase_index(PhysiCell_constants::live);
	cell_defaults.phenotype.cycle.data.transition_rate(cycle_start_index, cycle_end_index) = 0.0;

	// not motile 
	cell_defaults.phenotype.motility.is_motile = false;

	// add custom data here, if any 
	cell_defaults.custom_data.add_variable("time_of_death", "dimensionless", 0.0);
	cell_defaults.custom_data.add_variable("time_of_nucleation", "dimensionless", 0.0);

	return;
}

void setup_microenvironment(void)
{
	static int idxDeath = microenvironment.find_density_index("death_signal");

	default_microenvironment_options.simulate_2D = true;

	// make sure to override and go back to 2D 
	if (default_microenvironment_options.simulate_2D == false)
	{
		std::cout << "Warning: overriding XML config option and setting to 2D!" << std::endl;
		default_microenvironment_options.simulate_2D = true;
	}

	// override BioFVM setup with user parameters 
	microenvironment.diffusion_coefficients[idxDeath] = parameters.doubles("signal_diffusion_coefficient");

	initialize_microenvironment();

	return;
}

void setup_tissue(void)
{
	double length_x = microenvironment.mesh.bounding_box[3] -
		microenvironment.mesh.bounding_box[0];

	double length_y = microenvironment.mesh.bounding_box[4] -
		microenvironment.mesh.bounding_box[1];

	Cell* pC;

	std::string csv_file = parameters.strings("csv_file");
	std::string csv_file2 = parameters.strings("csv_file2");
	std::cout << "-------- reading " << csv_file << " and " << csv_file2 << std::endl;
	/*
	Cell index,Cell X,Cell Y,Time of Nucleation,Time of death
	1,2.25,6,,3
	2,3.3,4.25,1,1
	...
	*/
	std::ifstream infile(csv_file);
	std::ifstream infile2(csv_file2);
	static int idxDeath = microenvironment.find_density_index("death_signal");
	int idx, idx2;
	double x, y, t_nuc, t_death;
	char sep, sep2;
	std::string neighbors;
	int neiIdx;

	// Read a title header
	std::string line;
	std::string line2;
	std::getline(infile, line);
	std::getline(infile2, line2);
	std::cout << "-------- csv header: \n" << "id" << line << std::endl;
	static double radius = parameters.doubles("cell_radius");
	double volume = 4. / 3 * PhysiCell_constants::pi * radius * radius * radius;
	while ((infile >> idx >> sep >> x >> sep >> y >> sep >> t_nuc >> sep >> t_death) && (sep == ',') && std::getline(infile2, line2))
	{
		/*
		Cell index,Cell X,Cell Y,Time of Nucleation,Time of death
		1,2.25,6,,3
		2,3.3,4.25,1,1
		...
		*/
		std::cout << "cell " << idx << ":  x,y= " << x << "," << y;
		std::cout << ", t_nuc= " << t_nuc << ", t_death = " << t_death << std::endl;

		pC = create_cell();
		pC->index = idx;
		pC->assign_position(x, y, 0.0);
		pC->set_total_volume(volume);

		pC->phenotype.secretion.set_all_secretion_to_zero();
		pC->custom_data["time_of_death"] = 0.0;
		pC->custom_data["time_of_nucleation"] = t_nuc;
		pC->phenotype.molecular.internalized_total_substrates[idxDeath] = 0;

		//***just for test***
		if (idx >=90 && idx <= 100)
			pC->phenotype.molecular.internalized_total_substrates[idxDeath] = 1;
		
		//parse neighbors
		std::istringstream ss(line2);
		std::cout << line2 << std::endl;
		neiIdx = -1;
		while (ss) {
			std::string neighbor;
			if (!(std::getline(ss, neighbor, ',')) || neighbor.empty())
				break;
		
			try {
				if (neiIdx == -1) {
					idx2 = stoi(neighbor);
					std::cout << "Neighbors for cell  " << idx2 << std::endl;
					neiIdx++;
					std::getline(ss, neighbor, '"');  //for experiment real data only, to ignore "
					continue;
				}
				else {
					(*all_cells)[idx2]->state.neighbors.push_back((*all_cells)[stoi(neighbor)]);
					std::cout << "neighbor " << neiIdx << " = " << stoi(neighbor) << std::endl;
					neiIdx++;
				}
			}
			catch (const std::invalid_argument arg) {
				std::cout << "error in parsing the neighbors" << std::endl;
				arg.what();
			}
		}
	}

	return;
}

void death_function(Cell* pCell, Phenotype& phenotype, double dt)
{
	static int idxDeath = microenvironment.find_density_index("death_signal");
	double signal = phenotype.molecular.internalized_total_substrates[idxDeath];
	double nuc_time = pCell->custom_data["time_of_nucleation"];
	std::cout << "cell "<< pCell->index << "-------death_function:  i_t_s= " << signal << std::endl;
	std::cout << PhysiCell_globals.current_time << std::endl;

	if (signal >= parameters.doubles("death_threshold") || PhysiCell_globals.current_time >= nuc_time)
	{
		std::cout << "\t\tdie!" << std::endl;
		pCell->release_internalized_substrates();
		pCell->lyse_cell(); // start_death( apoptosis_model_index )
		return;
	}
	else if (signal >= parameters.doubles("min_death_count"))
	{
		double new_death = parameters.doubles("death_rate");
		new_death *= dt;
		std::cout << "cell " << pCell->index << "-------death_function:  increasing death= " << new_death << std::endl;
		phenotype.molecular.internalized_total_substrates[idxDeath] += new_death;
	}
		//absorb particles
		//double to_absorb = get_particles_from_surounding(pCell); //to complete function, checks the amount of particles in the cell's environment
		//double to_absorb = 0.1 * signal;
		//phenotype.molecular.internalized_total_substrates[idxDeath] += to_absorb;
		//remove to_absorb from environment
	
	return;

	// check for contact with a cell
	// Cell* pTestCell = NULL;
	// std::vector<Cell*> neighbors = get_possible_neighbors(pCell);

	// for( int n=0; n < neighbors.size() ; n++ )
	// {
	// 	pTestCell = neighbors[n];
	// 	// if it is not me and not a macrophage
	// 	if( pTestCell != pCell && pTestCell->type != macrophage.type )
	// 	{
	// 		// calculate distance to the cell
	// 		std::vector<double> displacement = pTestCell->position;
	// 		displacement -= pCell->position;
	// 		double distance = norm( displacement );

	// 		double max_distance = pCell->phenotype.geometry.radius +
	// 			pTestCell->phenotype.geometry.radius;
	// 		max_distance *= 1.1;

	// 		// if it is not a macrophage, test for viral load
	// 		// if high viral load, eat it.
	// 		if( pTestCell->phenotype.molecular.internalized_total_substrates[nSignal]
	// 			> parameters.doubles("min_virion_detection_threshold") &&
	// 			distance < max_distance )
	// 		{
	// 			std::cout << "\t\tnom nom nom" << std::endl;
	// 			pCell->ingest_cell( pTestCell );
	// 		}
	// 	}
	// }

}

std::vector<std::string> my_coloring_function(Cell* pCell)
{
	static int nSignal = microenvironment.find_density_index("death_signal");

	// start with flow cytometry coloring
	std::vector<std::string> output = false_cell_coloring_cytometry(pCell);

	std::cout << "--- my_coloring_function i_t_s " << pCell->phenotype.molecular.internalized_total_substrates[nSignal] << std::endl;

	if (pCell->phenotype.death.dead == false && pCell->type == 1)
	{
		output[0] = "black";
		output[2] = "black";
	}

	return output;
}


std::vector<std::string> death_coloring_function(Cell* pCell)
{
	// start with flow cytometry coloring 
	std::vector<std::string> output = { "magenta" , "black" , "magenta", "black" };
	static int idxDeath = microenvironment.find_density_index("death_signal");

	// static double min_virus = parameters.doubles( "min_virion_count" );
	// static double max_virus = parameters.doubles( "burst_virion_count" ); 
	// static double denominator = max_virus - min_virus + 1e-15; 

	std::cout << "cell: " << pCell->index <<"--- death_coloring:  cell i_t_s " << pCell->phenotype.molecular.internalized_total_substrates[idxDeath] << std::endl;
	output[0] = "cyan";
	output[1] = "cyan";
	if (pCell->phenotype.molecular.internalized_total_substrates[idxDeath] == 1) {
		output[0] = "white";
		output[1] = "white";
	}
	else if (pCell->phenotype.molecular.internalized_total_substrates[idxDeath] == 2) {
		output[0] = "green";
		output[1] = "green";
	}
	else {
		output[0] = "red";
		output[1] = "red";
	}

	// dead cells 
	// if( pCell->phenotype.death.dead == true )
	// {
	// 	 output[0] = "red"; 
	// 	 output[2] = "darkred"; 
	// 	 return output; 
	// }
	return output;
}

/*
std::vector <std::pair <int, int>> get_all_neighbors( )
{
	std::vector <std::pair<int, int> > neighbors = {};
	std::string csv_file = parameters.strings( "csv_file2" );
	std::cout << "-------- reading neighbors " << csv_file << std::endl;
	std::ifstream infile(csv_file);

	int idx1, idx2;
	char sep;

	std::string line;
	std::getline(infile, line);
	std::cout << "-------- csv header: \n" << line << std::endl;
	while ((infile >> idx1 >> sep >> idx2 ) && (sep == ','))
	{
		std::cout << "cell " << idx1 <<" is a neighbor of cell "<< idx2 << std::endl ;
		std::pair<int,int>  neighborCells= std::make_pair(idx1,idx2);
		neighbors.push_back(neighborCells);
	}

	return neighbors;
}


std::vector <Cell*> get_cell_neighbors(Cell* cell)
{
	std::vector<Cell*> neighbors;
	std::string csv_file = parameters.strings( "csv_file2" );
	std::cout << "" << std::endl;
	std::cout << "-------- reading neighbors info from" << csv_file << std::endl;
	std::ifstream infile(csv_file);

	int idx1, idx2;
	char sep;
	int currentIdx;
	std::string line;
	std::getline(infile, line);
	while ((infile >> idx1 >> sep >> idx2 ) && (sep == ','))
	{
		currentIdx=-1;
		if (cell->index == idx1){
			currentIdx=idx2;
		}
		else if (cell->index == idx2){
			currentIdx=idx1;
		 }
		if (currentIdx!=-1){ //find actual cell from all_cells
			for( unsigned int i=0; i < (*all_cells).size(); i++ )
				{
					Cell* pC = (*all_cells)[i];
					if (currentIdx == pC->index ){
						neighbors.push_back(pC);
						break;
					}
				 }
		}
	}

	return neighbors;
}
*/
std::vector<double> integrate_total_substrates(void)
{
	// start with 0 vector
	std::vector<double> out(microenvironment.number_of_densities(), 0.0);

	// integrate extracellular substrates
	for (unsigned int n = 0; n < microenvironment.number_of_voxels(); n++)
	{
		// out = out + microenvironment(n) * dV(n)
		axpy(&out, microenvironment.mesh.voxels[n].volume, microenvironment(n));
	}

	// inte
	for (unsigned int n = 0; n < (*all_cells).size(); n++)
	{
		Cell* pC = (*all_cells)[n];
		out += pC->phenotype.molecular.internalized_total_substrates;
	}

	return out;
}