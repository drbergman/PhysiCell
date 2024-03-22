#include "./extracellular_matrix.h"


extern ECM ecm;

double sign_function (double number)
{
	// double sign = 0.0
	if (number<0)
	{ return -1.0;}
	return 1.0;
}

void custom_update_motility_vector(Cell* pCell, Phenotype& phenotype, double dt_ )
{
	// Modified version of standard update_motility_vector function in PhysiCell. Search "update_motility_vector" to find original
	// Takes into account ECM anisotropy for anisotropy linked chemotactic migration
	// Also, outputs/modifies migration_bias_direction versus motility_vector (as in the original function)

	if( phenotype.motility.is_motile == false )
	{
		phenotype.motility.motility_vector.assign( 3, 0.0 ); 
		return; 
	}

	static bool is_first_time = true;
	is_first_time &= PhysiCell_globals.current_time == 0; // is_first_time starts true; switches to false the first time the current time is !=0 here; then stays false forever

	if (is_first_time || UniformRandom() < dt_ / phenotype.motility.persistence_time)
	{
		std::vector<double> randvec(3,0.0);
		if( phenotype.motility.restrict_to_2D == true )
		{ randvec = UniformOnUnitCircle(); }
		else
		{ randvec = UniformOnUnitSphere(); }
		
		// if the update_bias_vector function is set, use it  
		if( pCell->functions.update_migration_bias )
		{
			pCell->functions.update_migration_bias( pCell, phenotype, dt_ ); 
		}
		phenotype.motility.migration_bias_direction *= phenotype.motility.migration_bias; // motility = bias*bias_vector
		double one_minus_bias = 1.0 - phenotype.motility.migration_bias;
		axpy(&(phenotype.motility.migration_bias_direction), one_minus_bias, randvec); // motility = (1-bias)*randvec + bias*bias_vector
		normalize(&(phenotype.motility.migration_bias_direction));
	}
	return;
}

void ecm_to_cell_interactions( Cell* pCell, Phenotype& phenotype, double dt )
{
	static bool use_simple_fn = parameters.ints("ecm_model_id") == 1;
	if (use_simple_fn == true)
	{
		// Get cell level values
		static int ecm_sensitivity_index = pCell->custom_data.find_variable_index("ecm_sensitivity");
		if (ecm_sensitivity_index < 0)
		{
			std::cout << "        static int ecm_sensitivity_index = " << ecm_sensitivity_index << std::endl;
			std::exit(-1);
		}

		static int min_ecm_mot_den_index = pCell->custom_data.find_variable_index("min_ecm_motility_density");
		if (min_ecm_mot_den_index < 0)
		{
			std::cout << "        static int min_ecm_mot_den_index = " << min_ecm_mot_den_index << std::endl;
			std::exit(-1);
		}
		static int max_ecm_mot_den_index = pCell->custom_data.find_variable_index("max_ecm_motility_density");
		if (max_ecm_mot_den_index < 0)
		{
			std::cout << "        static int max_ecm_mot_den_index = " << max_ecm_mot_den_index << std::endl;
			std::exit(-1);
		}

		static int ideal_ecm_mot_den_index = pCell->custom_data.find_variable_index("ideal_ecm_motility_density");
		if (ideal_ecm_mot_den_index < 0)
		{
			std::cout << "        static int ideal_ecm_mot_den_index = " << ideal_ecm_mot_den_index << std::endl;
			std::exit(-1);
		}

		// sample ECM
		int nearest_ecm_voxel_index = ecm.ecm_mesh.nearest_voxel_index(pCell->position);
		double ecm_density = ecm.ecm_voxels[nearest_ecm_voxel_index].density;

		if (ecm_density <= pCell->custom_data[min_ecm_mot_den_index] || ecm_density >= pCell->custom_data[max_ecm_mot_den_index])
		{
			phenotype.motility.motility_vector = {0.0, 0.0, 0.0};
			return;
		}

		// Function calculates takes chemotaxis direction and combines with fiber following to produce combined random biased motility and
		// fiber following vector. Then calculates agent speed based on ECM density.

		std::vector<double> f = ecm.ecm_voxels[nearest_ecm_voxel_index].ecm_fiber_alignment;
		normalize(&f);

		// ******************************************************************************************************
		// ****** Make linear combination of random biased migration (motility_vector) and ECM orientation following. ********
		// ******************************************************************************************************

		std::vector<double> d_motility = phenotype.motility.migration_bias_direction;
		// normalize( &d_motility ); // not necessary since it is normalized at the end of custom_update_motility_vector

		// to determine direction along f, find part of d_choice that is perpendicular to f;
		std::vector<double> d_perp = d_motility - dot_product(d_motility, f) * f;

		normalize(&d_perp);

		// find constants to span d_choice with d_perp and f
		// double c_1 = dot_product(d_motility, d_perp);
		// double c_2 = dot_product(d_motility, f);

		// calculate bias away from directed motitility - combination of sensitity to ECM and anisotropy
		double gamma = pCell->custom_data[ecm_sensitivity_index] * ecm.ecm_voxels[nearest_ecm_voxel_index].anisotropy; // at low values, directed motility vector is recoved. At high values, fiber direction vector is recovered.

		phenotype.motility.motility_vector = (1.0 - gamma) * dot_product(d_motility, d_perp) * d_perp + dot_product(d_motility, f) * f;

		double directed_speed; // use this to set length of motility vector rather than phenotype.motility.migration_bias so that phenotype.motility.migration_bias can be saved between steps based on rules, etc. This will mean that the saved value of migration_speed will not account for the ecm effects

		if (parameters.bools("normalize_ecm_influenced_motility_vector") == true)
		{
			// if the vector is to be normalized, we, by definition, already know the magnitude will be 1.0
			directed_speed = 1.0;
		}
		else
		{
			directed_speed = norm(phenotype.motility.motility_vector);
		}
		normalize(&(phenotype.motility.motility_vector));

		/****************************************END new migration direction update****************************************/

		/*********************************************Begin speed update***************************************************/
		directed_speed *= pCell->phenotype.motility.migration_speed;
		double temp;
		// At this point, we have checked above that ecm_density \in [min,max] ecm motility density for this cell above
		if (ecm_density <= pCell->custom_data[ideal_ecm_mot_den_index])
		{
			temp = pCell->custom_data[min_ecm_mot_den_index];
		}
		else
		{
			temp = pCell->custom_data[max_ecm_mot_den_index];
		}
		directed_speed *= (ecm_density - temp) / (pCell->custom_data[ideal_ecm_mot_den_index] - temp);

		phenotype.motility.motility_vector *= directed_speed;
		/*********************************************END speed update***************************************************/
		return;
	}

	// else using more complicated function
	// Get cell level values
	static int ecm_sensitivity_density_ec50_at_anisotropy_0_index = pCell->custom_data.find_variable_index("ecm_sensitivity_density_ec50_at_anisotropy_0");
	static int ecm_sensitivity_density_ec50_at_anisotropy_1_index = pCell->custom_data.find_variable_index("ecm_sensitivity_density_ec50_at_anisotropy_1");
	if (pCell->custom_data[ecm_sensitivity_density_ec50_at_anisotropy_0_index] < pCell->custom_data[ecm_sensitivity_density_ec50_at_anisotropy_1_index] 
		|| pCell->custom_data[ecm_sensitivity_density_ec50_at_anisotropy_1_index] <= 0)
	{
		std::cout << "ERROR: The constraint 0 <= ecm_sensitivity_density_ec50_at_anisotropy_1 <= ecm_sensitivity_density_ec50_at_anisotropy_0 <= 1" << std::endl
				  << "\thas been violated." << std::endl;
		exit(-1);
	}
	else if (pCell->custom_data[ecm_sensitivity_density_ec50_at_anisotropy_0_index] > 1.0)
	{
		std::cout << "WARNING: The ecm_sensitivity_density_ec50_at_anisotropy_0 is greater than 1." << std::endl
				  << "\tDensity is expected to be on [0,1], so an EC50 beyond this range might not be what you want." << std::endl;
	}
	static int ecm_sensitivity_anisotropy_hill_coefficient_index = pCell->custom_data.find_variable_index("ecm_sensitivity_anisotropy_hill_coefficient");
	static int ecm_speed_increase_anisotropy_slope_index = pCell->custom_data.find_variable_index("ecm_speed_increase_anisotropy_slope");
	static int ecm_speed_increase_anisotropy_intercept_index = pCell->custom_data.find_variable_index("ecm_speed_increase_anisotropy_intercept");
	static int ecm_speed_increase_density_ec50_index = pCell->custom_data.find_variable_index("ecm_speed_increase_density_ec50");
	static int ecm_speed_increase_density_hill_coefficient_index = pCell->custom_data.find_variable_index("ecm_speed_increase_density_hill_coefficient");
	static int ecm_speed_decrease_anisotropy_slope_index = pCell->custom_data.find_variable_index("ecm_speed_decrease_anisotropy_slope");
	static int ecm_speed_decrease_anisotropy_intercept_index = pCell->custom_data.find_variable_index("ecm_speed_decrease_anisotropy_intercept");
	static int ecm_speed_decrease_density_hill_coefficient_index = pCell->custom_data.find_variable_index("ecm_speed_decrease_density_hill_coefficient");
	
	static int ecm_speed_factor_at_anisotropy_1_index = pCell->custom_data.find_variable_index("ecm_speed_factor_at_anisotropy_1");
	static int ecm_speed_factor_at_anisotropy_0_index = pCell->custom_data.find_variable_index("ecm_speed_factor_at_anisotropy_0");
	if (pCell->custom_data[ecm_speed_factor_at_anisotropy_1_index] < 0.0 || pCell->custom_data[ecm_speed_factor_at_anisotropy_1_index] > pCell->custom_data[ecm_speed_factor_at_anisotropy_0_index] || pCell->custom_data[ecm_speed_factor_at_anisotropy_0_index] > 1.0)
	{
		std::cout << " ERROR: The constraint 0 <= ecm_speed_factor_at_anisotropy_1 <= ecm_speed_factor_at_anisotropy_0 <= 1" << std::endl
				  << "\thas been violated." << std::endl;
		exit(-1);
	}

	static int ecm_min_speed_complex_ecm_model_index = pCell->custom_data.find_variable_index("ecm_min_speed_complex_ecm_model");
	// sample ECM
	int nearest_ecm_voxel_index = ecm.ecm_mesh.nearest_voxel_index(pCell->position);
	double ecm_density = ecm.ecm_voxels[nearest_ecm_voxel_index].density;

	// Function calculates takes chemotaxis direction and combines with fiber following to produce combined random biased motility and
	// fiber following vector. Then calculates agent speed based on ECM density.

	std::vector<double> f = ecm.ecm_voxels[nearest_ecm_voxel_index].ecm_fiber_alignment;
	normalize(&f);

	// ******************************************************************************************************
	// ****** Make linear combination of random biased migration (motility_vector) and ECM orientation following. ********
	// ******************************************************************************************************

	std::vector<double> d_motility = phenotype.motility.migration_bias_direction;
	// normalize( &d_motility ); // not necessary since it is normalized at the end of custom_update_motility_vector

	// to determine direction along f, find part of d_choice that is perpendicular to f;
	double cos_angle = dot_product(d_motility, f);
	std::vector<double> d_perp = d_motility - cos_angle * f;

	normalize(&d_perp);

	// find constants to span d_choice with d_perp and f
	// double c_1 = dot_product(d_motility, d_perp);
	// double c_2 = dot_product(d_motility, f);

	// calculate bias away from directed motitility - combination of sensitity to ECM and anisotropy
	double a = ecm.ecm_voxels[nearest_ecm_voxel_index].anisotropy;
	double ec50 = (1 - a) * pCell->custom_data[ecm_sensitivity_density_ec50_at_anisotropy_0_index] + a * pCell->custom_data[ecm_sensitivity_density_ec50_at_anisotropy_1_index];
	double n = pCell->custom_data[ecm_sensitivity_anisotropy_hill_coefficient_index];
	double sensitivity = a / (1 + pow(ecm_density/ec50,-n));

	phenotype.motility.motility_vector = (1.0 - sensitivity) * dot_product(d_motility, d_perp) * d_perp + dot_product(d_motility, f) * f;
	normalize(&(phenotype.motility.motility_vector));

	/****************************************END new migration direction update****************************************/

	/*********************************************Begin speed update***************************************************/
	// double directed_speed = pCell->phenotype.motility.migration_speed; // s_0
	double directed_speed = pCell->custom_data[ecm_min_speed_complex_ecm_model_index]; // s_0
	double ds = pCell->custom_data[ecm_speed_increase_anisotropy_slope_index] * a + pCell->custom_data[ecm_speed_increase_anisotropy_intercept_index];
	directed_speed += ds / (1 + pow(ecm_density / pCell->custom_data[ecm_speed_increase_density_ec50_index], -pCell->custom_data[ecm_speed_increase_density_hill_coefficient_index]));
	double gamma_2 = pCell->custom_data[ecm_speed_decrease_anisotropy_slope_index] * a + pCell->custom_data[ecm_speed_decrease_anisotropy_intercept_index];
	directed_speed /= 1 + pow(ecm_density / gamma_2, pCell->custom_data[ecm_speed_decrease_density_hill_coefficient_index]);
	// double f0 = pCell->custom_data[ecm_speed_factor_at_anisotropy_0_index] - (pCell->custom_data[ecm_speed_factor_at_anisotropy_0_index] - pCell->custom_data[ecm_speed_factor_at_anisotropy_1_index])*a;
	double f0 = pCell->custom_data[ecm_speed_factor_at_anisotropy_0_index] * (1 - a) + pCell->custom_data[ecm_speed_factor_at_anisotropy_1_index] * a;
	directed_speed *= f0 + (1 - f0) * fabs(cos_angle);

	phenotype.motility.motility_vector *= directed_speed;
	/*********************************************END speed update***************************************************/
	return;
}

void ecm_remodeling_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	// Find attributes needed for updating ECM

	int nearest_ecm_voxel_index = ecm.ecm_mesh.nearest_voxel_index( pCell->position );   
	static int Cell_ecm_production_rate_index = pCell->custom_data.find_variable_index( "ecm_production_rate");
	if (Cell_ecm_production_rate_index < 0) 
    {
        std::cout << "        static int Cell_ecm_production_rate_index = " <<Cell_ecm_production_rate_index << std::endl;
        std::exit(-1);   
    }
	static int cell_ecm_degradation_rate_index = pCell->custom_data.find_variable_index( "ecm_degradation_rate");
	if (cell_ecm_degradation_rate_index < 0) 
    {
        std::cout << "        static int cell_ecm_degradation_rate_index = " <<cell_ecm_degradation_rate_index << std::endl;
        std::exit(-1);   
    }
	static int cell_fiber_realignment_rate_index = pCell->custom_data.find_variable_index( "fiber_realignment_rate");
	if (cell_fiber_realignment_rate_index < 0) 
    {
        std::cout << "        static int cell_fiber_realignment_rate_index = " <<cell_fiber_realignment_rate_index << std::endl;
        std::exit(-1);
	}

	// Cell-ECM density interaction

    double ecm_density = ecm.ecm_voxels[nearest_ecm_voxel_index].density;
    double r = pCell->custom_data[Cell_ecm_production_rate_index];
    double d = pCell->custom_data[cell_ecm_degradation_rate_index];

	// ecm.ecm_voxels[nearest_ecm_voxel_index].density += dt * (r - d * ecm_density);
	if (d > 0)
	{
		double density_saturation = r / d;
		ecm.ecm_voxels[nearest_ecm_voxel_index].density = density_saturation + (ecm_density - density_saturation) * exp(-d * dt);
	}
	else
	{
		ecm.ecm_voxels[nearest_ecm_voxel_index].density += r*dt;
	}
	// End Cell-ECM density interaction

	// Cell-ECM Fiber realignment and anisotropy remodeling - continous then instantaneous

	std::vector<double> ecm_orientation = ecm.ecm_voxels[nearest_ecm_voxel_index].ecm_fiber_alignment;
	// double anisotropy = ecm.ecm_voxels[nearest_ecm_voxel_index].anisotropy;
	std::vector<double> norm_cell_velocity = pCell->velocity;
	double speed = norm(norm_cell_velocity);
	if (speed <= 1e-16)
	{
		norm_cell_velocity = {0.0, 0.0, 0.0};
	}
	else
	{
		for (int i = 0; i < 3; i++)
		{
			norm_cell_velocity[i] /= speed;
		}
		if (dot_product(norm_cell_velocity,ecm_orientation) < 0)
		for (int i = 0; i < 3; i++)
		{
			norm_cell_velocity[i] *= -1.0;
		}
	}
	double r_0 = pCell->custom_data[cell_fiber_realignment_rate_index] * speed;

	// normalize(&norm_cell_velocity);

	if (r_0 != 0 || r != 0)
	{ // only bother making these changes if the cell is actually changing the ECM structure
		for (int i = 0; i < 2; i++)
		{
			ecm.ecm_voxels[nearest_ecm_voxel_index].ecm_fiber_alignment[i] = norm_cell_velocity[i] - (norm_cell_velocity[i] - ecm_orientation[i]) * exp(-(r_0 + r / ecm_density) * dt);
		}
		ecm.ecm_voxels[nearest_ecm_voxel_index].anisotropy = norm(ecm.ecm_voxels[nearest_ecm_voxel_index].ecm_fiber_alignment);
	}
	// END Cell-ECM Anisotropy Modification
}

void custom_update_cell_velocity( Cell* pCell, Phenotype& phenotype, double dt)
{
	// Replaces the standard_update_cell_velocity in fiber following/senstive agents. 
	// Assign this to update_cell_velocity function pointer to get fiber following and density based speed changes
	

	if( pCell->functions.add_cell_basement_membrane_interactions )
	{
		pCell->functions.add_cell_basement_membrane_interactions(pCell, phenotype,dt);
	}
	
	pCell->state.simple_pressure = 0.0; 
	pCell->state.neighbors.clear(); 
	
	//First check the neighbors in my current voxel
	std::vector<Cell*>::iterator neighbor;
	std::vector<Cell*>::iterator end = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].end();
	for(neighbor = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].begin(); neighbor != end; ++neighbor)
	{
		pCell->add_potentials(*neighbor);
	}
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
		{
			pCell->add_potentials(*neighbor);
		}
	}

	// non-standard motility update - part of ECM package
	custom_update_motility_vector(pCell, phenotype, dt);

	// ECM following and speed update
	ecm_to_cell_interactions(pCell, phenotype, dt); 

	// standard update cell velocity - after this update proceeds as "conventional" PhysiCell
	pCell->velocity += phenotype.motility.motility_vector;

	ecm_remodeling_function(pCell, phenotype, dt);

	return; 
}

void create_default_ecm_compatible_agent( void )
{
	initialize_default_cell_definition(); 
	// cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.update_velocity = custom_update_cell_velocity;

	// cell_defaults.functions.custom_cell_rule = ecm_remodeling_function; 

	// 	cell_defaults.functions.cycle_model = Ki67_advanced; 
	
	// cell_defaults.functions.volume_update_function = standard_volume_update_function;
	// cell_defaults.functions.update_migration_bias = NULL; 
	
	// cell_defaults.functions.update_phenotype = update_cell_and_death_parameters_O2_based; // NULL; 
	// cell_defaults.functions.custom_cell_rule = NULL; 
	
	// cell_defaults.functions.update_velocity = standard_update_cell_velocity;
	// cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	// cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
}

void SVG_plot_custom( std::string filename , Microenvironment& M, double z_slice , double time, std::vector<std::string> (*cell_coloring_function)(Cell*), std::string line_pattern )
{
	double X_lower = M.mesh.bounding_box[0];
	double X_upper = M.mesh.bounding_box[3];
 
	double Y_lower = M.mesh.bounding_box[1]; 
	double Y_upper = M.mesh.bounding_box[4]; 

	double plot_width = X_upper - X_lower; 
	double plot_height = Y_upper - Y_lower; 

	double font_size = 0.025 * plot_height; // PhysiCell_SVG_options.font_size; 
	double top_margin = font_size*(.2+1+.2+.9+.5 ); 

	// open the file, write a basic "header"
	std::ofstream os( filename , std::ios::out );
	if( os.fail() )
	{ 
		std::cout << std::endl << "Error: Failed to open " << filename << " for SVG writing." << std::endl << std::endl; 

		std::cout << std::endl << "Error: We're not writing data like we expect. " << std::endl
		<< "Check to make sure your save directory exists. " << std::endl << std::endl
		<< "I'm going to exit with a crash code of -1 now until " << std::endl 
		<< "you fix your directory. Sorry!" << std::endl << std::endl; 
		exit(-1); 
	} 
	
	Write_SVG_start( os, plot_width , plot_height + top_margin );

	// draw the background 
	Write_SVG_rect( os , 0 , 0 , plot_width, plot_height + top_margin , 0.002 * plot_height , "white", "white" );
	
	// bool Write_SVG_circle( std::ostream& os, double center_x, double center_y, double radius, double stroke_size, 
    //                    std::string stroke_color , std::string fill_color )

	// write the simulation time to the top of the plot
 
	char* szString; 
	szString = new char [1024]; 
 
	int total_cell_count = all_cells->size(); 
 
	double temp_time = time; 

	std::string time_label = formatted_minutes_to_DDHHMM( temp_time ); 
 
	sprintf( szString , "Current time: %s, z = %3.2f %s", time_label.c_str(), 
		z_slice , PhysiCell_SVG_options.simulation_space_units.c_str() ); 
	Write_SVG_text( os, szString, font_size*0.5,  font_size*(.2+1), 
		font_size, PhysiCell_SVG_options.font_color.c_str() , PhysiCell_SVG_options.font.c_str() );
	sprintf( szString , "%u agents" , total_cell_count ); 
	Write_SVG_text( os, szString, font_size*0.5,  font_size*(.2+1+.2+.9), 
		0.95*font_size, PhysiCell_SVG_options.font_color.c_str() , PhysiCell_SVG_options.font.c_str() );
	
	delete [] szString; 

	
	// add an outer "g" for coordinate transforms 
	
	os << " <g id=\"tissue\" " << std::endl 
	   << "    transform=\"translate(0," << plot_height+top_margin << ") scale(1,-1)\">" << std::endl; 
	   
	// prepare to do mesh-based plot (later)
	
	double dx_stroma = M.mesh.dx; 
	double dy_stroma = M.mesh.dy; 
	
	os << "  <g id=\"ECM\">" << std::endl; 
  
	int ratio = 1; 
	double voxel_size = dx_stroma / (double) ratio ; 
  
	double half_voxel_size = voxel_size / 2.0; 
	double normalizer = 78.539816339744831 / (voxel_size*voxel_size*voxel_size); 
 
 // color in the background ECM
/* 
 if( ECM.TellRows() > 0 )
 {
  // find the k corresponding to z_slice
  
  
  
  Vector position; 
  *position(2) = z_slice; 
  

  // 25*pi* 5 microns^2 * length (in source) / voxelsize^3
  
  for( int j=0; j < ratio*ECM.TellCols() ; j++ )
  {
   // *position(1) = *Y_environment(j); 
   *position(1) = *Y_environment(0) - dy_stroma/2.0 + j*voxel_size + half_voxel_size; 
   
   for( int i=0; i < ratio*ECM.TellRows() ; i++ )
   {
    // *position(0) = *X_environment(i); 
    *position(0) = *X_environment(0) - dx_stroma/2.0 + i*voxel_size + half_voxel_size; 
	
    double E = evaluate_Matrix3( ECM, X_environment , Y_environment, Z_environment , position );	
	double BV = normalizer * evaluate_Matrix3( OxygenSourceHD, X_environment , Y_environment, Z_environment , position );
	if( isnan( BV ) )
	{ BV = 0.0; }

	vector<string> Colors;
	Colors = hematoxylin_and_eosin_stroma_coloring( E , BV );
	Write_SVG_rect( os , *position(0)-half_voxel_size-X_lower , *position(1)-half_voxel_size+top_margin-Y_lower, 
	voxel_size , voxel_size , 1 , Colors[0], Colors[0] );
   
   }
  }
 
 }
*/
	os << "  </g>" << std::endl; 
 
	// Now draw vessels

	/*
	 std::vector<std::string> VesselColors = hematoxylin_and_eosin_stroma_coloring( 0,1 );

	 os << " <g id=\"BloodVessels\">" << endl; 
	 extern vector<BloodVesselSegment*> BloodVesselSegments; 
	 Vector Offset; 
	 *Offset(0) = X_lower; 
	 *Offset(1) = Y_lower-top_margin;
	*/
 

 
	// plot intersecting cells 
	os << "  <g id=\"cells\">" << std::endl; 
	for( int i=0 ; i < total_cell_count ; i++ )
	{
		Cell* pC = (*all_cells)[i]; // global_cell_list[i]; 
  
		static std::vector<std::string> Colors; 
		if( fabs( (pC->position)[2] - z_slice ) < pC->phenotype.geometry.radius )
		{
			double r = pC->phenotype.geometry.radius ; 
			double rn = pC->phenotype.geometry.nuclear_radius ; 
			double z = fabs( (pC->position)[2] - z_slice) ; 
   
			Colors = cell_coloring_function( pC ); 

			os << "   <g id=\"cell" << pC->ID << "\">" << std::endl; 
  
			// figure out how much of the cell intersects with z = 0 
   
			double plot_radius = sqrt( r*r - z*z ); 

			Write_SVG_circle( os, (pC->position)[0]-X_lower, (pC->position)[1]-Y_lower, 
				plot_radius , 0.5, Colors[1], Colors[0] ); 

			// plot the nucleus if it, too intersects z = 0;
			if( fabs(z) < rn && PhysiCell_SVG_options.plot_nuclei == true )
			{   
				plot_radius = sqrt( rn*rn - z*z ); 
			 	Write_SVG_circle( os, (pC->position)[0]-X_lower, (pC->position)[1]-Y_lower, 
					plot_radius, 0.5, Colors[3],Colors[2]); 
			}					  
			os << "   </g>" << std::endl;
		}
	}
	os << "  </g>" << std::endl; 

	// Plot guidelines the circle guides
	
	if (line_pattern == "concentric circles")
	{
		for (int i=0; i<plot_width/(2*parameters.doubles("ecm_dx")); i++)
		{
			double radius = plot_width/2 - (i *parameters.doubles("ecm_dx"));
			// std::cout<<"Index "<<i<<std::endl;
			// std::cout<<"Radius "<<radius<<std::endl;
			Write_SVG_circle( os, plot_width/2, plot_height/2, radius, 0.5, "black" , "none" );

		}
	}

	else if (line_pattern == "vertical lines")
	{
		for (int i=0; i<plot_width/parameters.doubles("ecm_dx"); i++)
		{
			double x_line_position = parameters.doubles("ecm_dx")*i;
			// std::cout<<"Index "<<i<<std::endl;
			// std::cout<<"X position "<<x_line_position<<std::endl;			
			Write_SVG_line(os, x_line_position, 0, x_line_position, plot_height, 0.5, "black");
			// bool Write_SVG_line( std::ostream& os , double start_x, double start_y, double end_x , double end_y, double thickness, 
            //         std::string stroke_color )
		}
	}

	else if (line_pattern == "horizontal lines")
	{
		for (int i=0; i<plot_height/parameters.doubles("ecm_dy"); i++)
		{
			double y_line_position = parameters.doubles("ecm_dy")*i;
			// std::cout<<"Index "<<i<<std::endl;
			// std::cout<<"Y position "<<y_line_position<<std::endl;			
			Write_SVG_line(os, 0, y_line_position, plot_width, y_line_position, 0.5, "black");
			// bool Write_SVG_line( std::ostream& os , double start_x, double start_y, double end_x , double end_y, double thickness, 
            //         std::string stroke_color )
		}
	}

	else if (line_pattern == "none") {} // Don't make lines!!!

	else if (line_pattern != "none" || "horizontal lines" || "vertical lines" || "concentric circles")
	{
		std::cout<<"Use of this custom SVG output function requires specifying \"none\" or a line pattern" <<std::endl;
		std::cout<<"Halting: see inputs to custom SVG function \"SVG_plot_custom\"" << std::endl;
		abort();
		return;

	}


	// plot intersecting BM points
	/* 
	 for( int i=0 ; i < BasementMembraneNodes.size() ; i++ )
	 {
		// vector<string> Colors = false_cell_coloring( pC ); 
		BasementMembraneNode* pBMN = BasementMembraneNodes[i]; 
		double thickness =0.1; 
		
		if( fabs( *(pBMN->Position)(2) - z_slice ) < thickness/2.0 ) 
		{
		 string bm_color ( "rgb(0,0,0)" );
		 double r = thickness/2.0; 
		 double z = fabs( *(pBMN->Position)(2) - z_slice) ; 

		 os << " <g id=\"BMN" << pBMN->ID << "\">" << std::endl; 
		 Write_SVG_circle( os,*(pBMN->Position)(0)-X_lower, *(pBMN->Position)(1)+top_margin-Y_lower, 10*thickness/2.0 , 0.5 , bm_color , bm_color ); 
		 os << " </g>" << std::endl;
		}
		// pC = pC->pNextCell;
	 }
	*/ 
	
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

	if( temp > 999 && std::strstr( bar_units.c_str() , PhysiCell_SVG_options.mu.c_str() )   )
	{
		temp /= 1000;
		bar_units = "mm";
	}
	// convert from mm to cm 
	if( temp > 9 && std::strcmp( bar_units.c_str() , "mm" ) == 0 )
	{
		temp /= 10; 
		bar_units = "cm";
	}
	
	szString = new char [1024];
	sprintf( szString , "%u %s" , (int) round( temp ) , bar_units.c_str() );
 
	Write_SVG_rect( os , plot_width - bar_margin - bar_width  , plot_height + top_margin - bar_margin - bar_height , 
		bar_width , bar_height , 0.002 * plot_height , "rgb(255,255,255)", "rgb(0,0,0)" );
	Write_SVG_text( os, szString , plot_width - bar_margin - bar_width + 0.25*font_size , 
		plot_height + top_margin - bar_margin - bar_height - 0.25*font_size , 
		font_size , PhysiCell_SVG_options.font_color.c_str() , PhysiCell_SVG_options.font.c_str() ); 
	
	delete [] szString; 

	// // plot runtime 
	// szString = new char [1024]; 
	// RUNTIME_TOC(); 
	// std::string formatted_stopwatch_value = format_stopwatch_value( runtime_stopwatch_value() );
	// Write_SVG_text( os, formatted_stopwatch_value.c_str() , bar_margin , top_margin + plot_height - bar_margin , 0.75 * font_size , 
	// 	PhysiCell_SVG_options.font_color.c_str() , PhysiCell_SVG_options.font.c_str() );
	// delete [] szString; 

	// draw a box around the plot window
	Write_SVG_rect( os , 0 , top_margin, plot_width, plot_height , 0.002 * plot_height , "rgb(0,0,0)", "none" );
	
	// close the svg tag, close the file
	Write_SVG_end( os ); 
	os.close();
 
	return; 
}

void write_ecm_Data_matlab( std::string filename )
{
    int number_of_data_entries = ecm.ecm_mesh.voxels.size();
    int size_of_each_datum = 8;
	
	// static int ecm_anisotropy_index = microenvironment.find_density_index( "ECM anisotropy" ); 
	// static int ecm_density_index = microenvironment.find_density_index( "ECM" ); 

    FILE* fp = write_matlab_header( size_of_each_datum, number_of_data_entries,  filename, "ecm_Data" );  // Note - the size of datum needs to correspond exaectly to the lines of output or there is an error upon importing.

    for( int i=0; i < number_of_data_entries ; i++ )
    {
	    fwrite( (char*) &( ecm.ecm_mesh.voxels[i].center[0] ) , sizeof(double) , 1 , fp ); // 1
        fwrite( (char*) &( ecm.ecm_mesh.voxels[i].center[1] ) , sizeof(double) , 1 , fp ); // 2
        fwrite( (char*) &( ecm.ecm_mesh.voxels[i].center[2] ) , sizeof(double) , 1 , fp ); //3
		fwrite( (char*) &( ecm.ecm_voxels[i].anisotropy), sizeof(double) , 1 , fp ); // 4
        fwrite( (char*) &( ecm.ecm_voxels[i].density), sizeof(double) , 1 , fp ); // 5
        fwrite( (char*) &( ecm.ecm_voxels[i].ecm_fiber_alignment[0]), sizeof(double) , 1 , fp ); // 6
        fwrite( (char*) &( ecm.ecm_voxels[i].ecm_fiber_alignment[1]), sizeof(double) , 1 , fp ); // 7
        fwrite( (char*) &( ecm.ecm_voxels[i].ecm_fiber_alignment[2]), sizeof(double) , 1 , fp ); // 8

		// This will only work if the diffusion and ECM meshes are the same size. Commenting out for actualrunning. To do a direct comparison, leave them in and change length. Will have to change the vizualization to get this from the regular BioFVM outputs.

		// fwrite( (char*) &( microenvironment.gradient_vector(i)[0][0]), sizeof(double) , 1 , fp ); // 9
		// fwrite( (char*) &( microenvironment.gradient_vector(i)[0][1]), sizeof(double) , 1 , fp ); // 10
		// fwrite( (char*) &( microenvironment.gradient_vector(i)[0][2]), sizeof(double) , 1 , fp ); // 11

    }
    fclose( fp );
    return;
}