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

std::ofstream deadfile;
extern std::vector<double> x_dead, y_dead;

// declare cell definitions here 
void create_cell_types(void)
{
	deadfile.open("dead_cells.dat");

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
	cell_defaults.custom_data.add_variable("original_index", "dimensionless", 0.0);

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
	while ((infile >> idx >> sep >> x >> sep >> y >> sep >> t_death >> sep >> t_nuc) && (sep == ',') && std::getline(infile2, line2))
	{
		std::cout << "cell " << idx << ":  x,y= " << x << "," << y;
		std::cout << ", t_death= " << t_death << ", t_nuc = " << t_nuc << std::endl;

		pC = create_cell();
		pC->index = idx;
		pC->assign_position(x, y, 0.0);
		pC->set_total_volume(volume);

		pC->phenotype.secretion.set_all_secretion_to_zero();
		pC->custom_data["time_of_death"] = t_death * 15;
		pC->custom_data["time_of_nucleation"] = t_nuc * 15;
		pC->custom_data["original_index"] = idx;
		pC->phenotype.molecular.internalized_total_substrates[idxDeath] = parameters.doubles("initial_internalized_substrate");

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
	double nuc_time = pCell->custom_data["time_of_nucleation"]; // real time of nucleation
	std::cout << "cell "<< pCell->index << "-------death_function:  i_t_s= " << signal << std::endl;
	std::cout << "time: " << PhysiCell_globals.current_time << std::endl;

	if (signal >= parameters.doubles("death_threshold") || PhysiCell_globals.current_time >= nuc_time)
	{
		//write dead cell info to file
		deadfile << PhysiCell_globals.current_time << " " << pCell->custom_data["time_of_death"] << " " << pCell->custom_data["original_index"] << " " << pCell->position[0] << "," << pCell->position[1] << std::endl;
		//std::cout << "\t\tdie!" << std::endl
		std::cout << "-------- save dead cell at " << pCell->position[0] << ',' << pCell->position[1] << std::endl;
		x_dead.push_back(pCell->position[0]);
		y_dead.push_back(pCell->position[1]);
		pCell->release_internalized_substrates();
		pCell->lyse_cell();                   // start_death( apoptosis_model_index )	
	}
	
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

	return output;
}


std::vector <Cell*> get_cell_neighbors(Cell* pCell)
{
	return pCell->state.neighbors;
}


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

//---------------------------
void my_SVG_plot(std::string filename, Microenvironment& M, double z_slice, double time, std::vector<std::string>(*cell_coloring_function)(Cell*))
{
	double X_lower = M.mesh.bounding_box[0];
	double X_upper = M.mesh.bounding_box[3];

	double Y_lower = M.mesh.bounding_box[1];
	double Y_upper = M.mesh.bounding_box[4];

	double plot_width = X_upper - X_lower;
	double plot_height = Y_upper - Y_lower;

	double font_size = 0.025 * plot_height; // PhysiCell_SVG_options.font_size; 
	double top_margin = font_size * (.2 + 1 + .2 + .9 + .5);

	// open the file, write a basic "header"
	std::ofstream os(filename, std::ios::out);
	if (os.fail())
	{
		std::cout << std::endl << "Error: Failed to open " << filename << " for SVG writing." << std::endl << std::endl;

		std::cout << std::endl << "Error: We're not writing data like we expect. " << std::endl
			<< "Check to make sure your save directory exists. " << std::endl << std::endl
			<< "I'm going to exit with a crash code of -1 now until " << std::endl
			<< "you fix your directory. Sorry!" << std::endl << std::endl;
		exit(-1);
	}

	Write_SVG_start(os, plot_width, plot_height + top_margin);

	// draw the background 
	Write_SVG_rect(os, 0, 0, plot_width, plot_height + top_margin, 0.002 * plot_height, "white", "white");

	// write the simulation time to the top of the plot

	char* szString;
	szString = new char[1024];

	int total_cell_count = all_cells->size();

	double temp_time = time;

	std::string time_label = formatted_minutes_to_DDHHMM(temp_time);

	sprintf(szString, "Current time: %s, z = %3.2f %s", time_label.c_str(),
		z_slice, PhysiCell_SVG_options.simulation_space_units.c_str());
	Write_SVG_text(os, szString, font_size * 0.5, font_size * (.2 + 1),
		font_size, PhysiCell_SVG_options.font_color.c_str(), PhysiCell_SVG_options.font.c_str());
	sprintf(szString, "%u agents", total_cell_count);
	Write_SVG_text(os, szString, font_size * 0.5, font_size * (.2 + 1 + .2 + .9),
		0.95 * font_size, PhysiCell_SVG_options.font_color.c_str(), PhysiCell_SVG_options.font.c_str());

	delete[] szString;


	// add an outer "g" for coordinate transforms 

	os << " <g id=\"tissue\" " << std::endl
		<< "    transform=\"translate(0," << plot_height + top_margin << ") scale(1,-1)\">" << std::endl;

	// prepare to do mesh-based plot (later)

	double dx_stroma = M.mesh.dx;
	double dy_stroma = M.mesh.dy;

	os << "  <g id=\"ECM\">" << std::endl;

	int ratio = 1;
	double voxel_size = dx_stroma / (double)ratio;

	double half_voxel_size = voxel_size / 2.0;
	double normalizer = 78.539816339744831 / (voxel_size * voxel_size * voxel_size);
	os << "  </g>" << std::endl;

	// plot intersecting cells 
	os << "  <g id=\"cells\">" << std::endl;
	for (int i = 0; i < total_cell_count; i++)
	{
		Cell* pC = (*all_cells)[i]; // global_cell_list[i]; 

		static std::vector<std::string> Colors;
		if (fabs((pC->position)[2] - z_slice) < pC->phenotype.geometry.radius)
		{
			double r = pC->phenotype.geometry.radius;
			double rn = pC->phenotype.geometry.nuclear_radius;
			double z = fabs((pC->position)[2] - z_slice);

			Colors = cell_coloring_function(pC);

			os << "   <g id=\"cell" << pC->ID << "\">" << std::endl;

			// figure out how much of the cell intersects with z = 0 

			double plot_radius = sqrt(r * r - z * z);

			Write_SVG_circle(os, (pC->position)[0] - X_lower, (pC->position)[1] - Y_lower,
				plot_radius, 0.5, Colors[1], Colors[0]);

			// plot the nucleus if it, too intersects z = 0;
			if (fabs(z) < rn && PhysiCell_SVG_options.plot_nuclei == true)
			{
				plot_radius = sqrt(rn * rn - z * z);
				Write_SVG_circle(os, (pC->position)[0] - X_lower, (pC->position)[1] - Y_lower,
					plot_radius, 0.5, Colors[3], Colors[2]);
			}
			os << "   </g>" << std::endl;
		}
	}

	// plot dead cells
/* e.g.,
<g id="dead001">
  <circle cx="22.018" cy="1352.51" r="20" stroke-width="0.5" stroke="red" fill="red"/>
</g>
*/
	double thickness = 0.1;
	std::string dead_color("rgb(100,100,100)");
	for (int idx = 0; idx < x_dead.size(); idx++)
	{
		os << " <g id=\"dead" << std::to_string(idx) << "\">" << std::endl;
		Write_SVG_circle(os, x_dead[idx] - X_lower, y_dead[idx] - Y_lower, 20, thickness, dead_color, dead_color);
		// os << "  <circle cx=\"" << std::to_string(x_dead[idx]) << "\"" <<  " cy=\"" << std::to_string(y_dead[idx]) << "\"  r=\"20\"" <<  
		// y_dead[idx], thickness, 0.5 , dead_color , dead_color ); 
		os << " </g>" << std::endl;
	}

	os << "  </g>" << std::endl;

	// end of the <g ID="tissue">
	os << " </g>" << std::endl;

	// draw a scale bar

	double bar_margin = 0.025 * plot_height;
	double bar_height = 0.01 * plot_height;
	double bar_width = PhysiCell_SVG_options.length_bar;
	double bar_stroke_width = 0.001 * plot_height;

	std::string bar_units = PhysiCell_SVG_options.simulation_space_units;
	// convert from micron to mm
	double temp = bar_width;

	if (temp > 999 && std::strstr(bar_units.c_str(), PhysiCell_SVG_options.mu.c_str()))
	{
		temp /= 1000;
		bar_units = "mm";
	}
	// convert from mm to cm 
	if (temp > 9 && std::strcmp(bar_units.c_str(), "mm") == 0)
	{
		temp /= 10;
		bar_units = "cm";
	}

	szString = new char[1024];
	sprintf(szString, "%u %s", (int)round(temp), bar_units.c_str());

	Write_SVG_rect(os, plot_width - bar_margin - bar_width, plot_height + top_margin - bar_margin - bar_height,
		bar_width, bar_height, 0.002 * plot_height, "rgb(255,255,255)", "rgb(0,0,0)");
	Write_SVG_text(os, szString, plot_width - bar_margin - bar_width + 0.25 * font_size,
		plot_height + top_margin - bar_margin - bar_height - 0.25 * font_size,
		font_size, PhysiCell_SVG_options.font_color.c_str(), PhysiCell_SVG_options.font.c_str());

	delete[] szString;

	// plot runtime 
	szString = new char[1024];
	RUNTIME_TOC();
	std::string formatted_stopwatch_value = format_stopwatch_value(runtime_stopwatch_value());
	Write_SVG_text(os, formatted_stopwatch_value.c_str(), bar_margin, top_margin + plot_height - bar_margin, 0.75 * font_size,
		PhysiCell_SVG_options.font_color.c_str(), PhysiCell_SVG_options.font.c_str());
	delete[] szString;

	// draw a box around the plot window
	Write_SVG_rect(os, 0, top_margin, plot_width, plot_height, 0.002 * plot_height, "rgb(0,0,0)", "none");

	// close the svg tag, close the file
	Write_SVG_end(os);
	os.close();

	return;
}