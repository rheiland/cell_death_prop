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

void create_cell_types( void )
{
	SeedRandom( parameters.ints("random_seed") ); // or specify a seed here 

	// housekeeping 
	initialize_default_cell_definition();

	// Name the default cell type 
	cell_defaults.type = 0; 
	cell_defaults.name = "my cell"; 
	cell_defaults.phenotype.mechanics.cell_cell_adhesion_strength = 0.0; 

	// set default cell cycle model 
	// cell_defaults.functions.cycle_model = flow_cytometry_separated_cycle_model; 

	// set default_cell_functions; 
	// cell_defaults.functions.update_phenotype = epithelial_function;

	// needed for a 2-D simulation: 
	
	cell_defaults.functions.set_orientation = up_orientation; 
	cell_defaults.phenotype.geometry.polarity = 1.0;
	cell_defaults.phenotype.motility.restrict_to_2D = true; 

	// make sure the defaults are self-consistent. 

	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	cell_defaults.phenotype.molecular.sync_to_microenvironment( &microenvironment );

	cell_defaults.functions.update_phenotype = death_function; 
	cell_defaults.phenotype.sync_to_functions( cell_defaults.functions ); 

	// set the rate terms in the default phenotype 

	// first find index for a few key variables. 
	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	int necrosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Necrosis" );
	// int virus_index = microenvironment.find_density_index( "virus" ); 

	// int G0G1_index = flow_cytometry_separated_cycle_model.find_phase_index( PhysiCell_constants::G0G1_phase );
	// int S_index = flow_cytometry_separated_cycle_model.find_phase_index( PhysiCell_constants::S_phase );

	// initially no death 
	cell_defaults.phenotype.death.rates[apoptosis_model_index] = 0.0; 
	cell_defaults.phenotype.death.rates[necrosis_model_index] = 0.0; 

	// initially no birth 
	// cell_defaults.phenotype.cycle.data.transition_rate(G0G1_index, S_index ) = 0.0 ; 

	// not motile 
	cell_defaults.phenotype.motility.is_motile = false; 

	// add custom data here, if any 
	cell_defaults.custom_data.add_variable( "time_of_death", "dimensionless", 0.0 );

	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
/* now this is in XML 
	default_microenvironment_options.X_range = {-1000, 1000}; 
	default_microenvironment_options.Y_range = {-1000, 1000}; 
	default_microenvironment_options.simulate_2D = true; 
*/
	// make sure to override and go back to 2D 
	if( default_microenvironment_options.simulate_2D == false )
	{
		std::cout << "Warning: overriding XML config option and setting to 2D!" << std::endl; 
		default_microenvironment_options.simulate_2D = true; 
	}
	
	// override BioFVM setup with user parameters 
	// int dummy_ID = microenvironment.find_density_index( "dummy" ); 
	int signal_ID = microenvironment.find_density_index( "death_signal" ); 
	// microenvironment.diffusion_coefficients[signal_ID] = parameters.doubles("viral_diffusion_coefficient"); 
	
	default_microenvironment_options.track_internalized_substrates_in_each_agent = true; 
	default_microenvironment_options.outer_Dirichlet_conditions = false;
	microenvironment.diffusion_coefficients[0] = 10; 	// no diffusion

	// initialize BioFVM 
	initialize_microenvironment(); 	

	// double factor = 1.0;
	// for( unsigned int n=0; n < microenvironment.number_of_voxels() ; n++ )
	// {
	// 	factor = 1.0 - ((n%20)/20.);
	// 	bc_vector[0] = factor*38.;
	// 	microenvironment(n) = bc_vector; 
	// }	

	// std::vector<double> bc_vector( 3 , 42  ); // 21% o2
	// bc_vector[1] = 1.0; 
	// bc_vector[2] = 0.0; 


	std::vector<double> bc_vector( 1 , 42  );   // 21% o2
	int n = microenvironment.number_of_voxels() / 2;
	bc_vector[0] = 500.0; 
	microenvironment(n+20) = bc_vector; 
	
	return; 
}

void setup_tissue( void )
{
	double length_x = microenvironment.mesh.bounding_box[3] - 
		microenvironment.mesh.bounding_box[0]; 
		
	double length_y = microenvironment.mesh.bounding_box[4] - 
		microenvironment.mesh.bounding_box[1]; 
		
	Cell* pC;
	
	std::string csv_file = parameters.strings( "csv_file" ); 
	std::cout << "-------- reading " << csv_file << std::endl;
/*
Cell index,Cell X,Cell Y,Time of Nucleation,Time of death
1,2.25,6,,3
2,3.3,4.25,1,1
...
*/
	std::ifstream infile(csv_file);

	static int nSignal = microenvironment.find_density_index( "death_signal" ); 
	int idx;
	double x,y, t_nuc,t_death;
	char sep;

	std::string line;
	std::getline(infile, line);
	std::cout << "-------- csv header: \n" << line << std::endl;
	// while (std::getline(infile, line) && (c == ',') )
	static double radius = parameters.doubles( "cell_radius" );
	double volume = 4./3 * PhysiCell_constants::pi * radius*radius*radius;
	while ((infile >> idx >> sep >> x >> sep >> y >> sep >> t_nuc >> sep >> t_death) && (sep == ','))
	{
/*
Cell index,Cell X,Cell Y,Time of Nucleation,Time of death
1,2.25,6,,3
2,3.3,4.25,1,1
...
*/
		std::cout << "cell " << idx <<":  x,y= "<< x << "," << y ;
		std::cout << ", t_nuc= " << t_nuc << ", t_death = " << t_death << std::endl;

		pC = create_cell(); 
		pC->assign_position( x,y, 0.0 );
		pC->set_total_volume(volume);

		pC->phenotype.secretion.set_all_secretion_to_zero(); 
		// pC->functions.update_phenotype = NULL; 
		pC->custom_data[ "time_of_death" ] = 0.0; 
		
		pC->phenotype.molecular.internalized_total_substrates[ nSignal ] = 1; 
	}

	return; 
}

void death_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int nSignal = microenvironment.find_density_index( "death_signal" ); 
	
	// digest virus particles inside me 
	static double implicit_Euler_constant = 
		(1.0 + dt * parameters.doubles("virus_digestion_rate") );  // todo: check for defined!

	// phenotype.molecular.internalized_total_substrates[nSignal] /= implicit_Euler_constant; 
	// phenotype.molecular.internalized_total_substrates[nSignal] /= implicit_Euler_constant; 
	phenotype.molecular.internalized_total_substrates[nSignal] += 1; 
	std::cout << "---------  int_tot_sub= " << phenotype.molecular.internalized_total_substrates[nSignal] << std::endl;
	
	if (phenotype.molecular.internalized_total_substrates[nSignal] > 2) 
	{
		std::cout << "\t\tdie!" << std::endl; 
		pCell->lyse_cell(); // start_death( apoptosis_model_index );
		pCell->functions.update_phenotype = NULL; 
		return; 
	}

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
	
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	static int nSignal = microenvironment.find_density_index( "death_signal" ); 
	
	// start with flow cytometry coloring 
	std::vector<std::string> output = false_cell_coloring_cytometry(pCell); 

	std::cout << "--- my_coloring_function i_t_s " <<  pCell->phenotype.molecular.internalized_total_substrates[nSignal] << std::endl;
		
	if( pCell->phenotype.death.dead == false && pCell->type == 1 )
	{
		 output[0] = "black"; 
		 output[2] = "black"; 
	}
	
	return output; 
}

std::vector<std::string> death_coloring_function( Cell* pCell )
{
	// start with flow cytometry coloring 
	
	std::vector<std::string> output = { "magenta" , "black" , "magenta", "black" }; 
	static int nSignal = microenvironment.find_density_index( "death_signal" ); 
	
	// static double min_virus = parameters.doubles( "min_virion_count" );
	// static double max_virus = parameters.doubles( "burst_virion_count" ); 
	// static double denominator = max_virus - min_virus + 1e-15; 

	std::cout << "--- cell i_t_s " <<  pCell->phenotype.molecular.internalized_total_substrates[nSignal] << std::endl;

	output[0] = "black"; 
	if( pCell->phenotype.molecular.internalized_total_substrates[nSignal] == 1 )
		 output[0] = "red"; 
	else if( pCell->phenotype.molecular.internalized_total_substrates[nSignal] == 2 )
		 output[0] = "green"; 
	else 
		 output[0] = "blue"; 
				
	// dead cells 
	// if( pCell->phenotype.death.dead == true )
	// {
	// 	 output[0] = "red"; 
	// 	 output[2] = "darkred"; 
	// 	 return output; 
	// }
	return output; 
}

std::vector<Cell*> get_possible_neighbors( Cell* pCell )
{
	std::vector<Cell*> neighbors = {}; 

	// First check the neighbors in my current voxel
	std::vector<Cell*>::iterator neighbor;
	std::vector<Cell*>::iterator end =
		pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].end();
	for( neighbor = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].begin(); neighbor != end; ++neighbor)
	{ neighbors.push_back( *neighbor ); }

	std::vector<int>::iterator neighbor_voxel_index;
	std::vector<int>::iterator neighbor_voxel_index_end = 
		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].end();

	for( neighbor_voxel_index = 
		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].begin();
		neighbor_voxel_index != neighbor_voxel_index_end; 
		++neighbor_voxel_index )
	{
		if(!is_neighbor_voxel(pCell, pCell->get_container()->underlying_mesh.voxels[pCell->get_current_mechanics_voxel_index()].center, pCell->get_container()->underlying_mesh.voxels[*neighbor_voxel_index].center, *neighbor_voxel_index))
			continue;
		end = pCell->get_container()->agent_grid[*neighbor_voxel_index].end();
		for(neighbor = pCell->get_container()->agent_grid[*neighbor_voxel_index].begin();neighbor != end; ++neighbor)
		{ neighbors.push_back( *neighbor ); }
	}
	
	return neighbors; 
}

std::vector<double> integrate_total_substrates( void )
{
	// start with 0 vector 
	std::vector<double> out( microenvironment.number_of_densities() , 0.0 ); 

	// integrate extracellular substrates 
	for( unsigned int n = 0; n < microenvironment.number_of_voxels() ; n++ )
	{
		// out = out + microenvironment(n) * dV(n) 
		axpy( &out , microenvironment.mesh.voxels[n].volume , microenvironment(n) ); 
	}

	// inte
	for( unsigned int n=0; n < (*all_cells).size(); n++ )
	{
		Cell* pC = (*all_cells)[n];
		out += pC->phenotype.molecular.internalized_total_substrates;
	}
	
	return out; 
}