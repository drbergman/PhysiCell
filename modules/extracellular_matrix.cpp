#include "./extracellular_matrix.h"

ECM ecm;

ECM_Cartesian_Mesh::ECM_Cartesian_Mesh()
{

	Cartesian_Mesh();
}

ECM_Voxel::ECM_Voxel()
{
	anisotropy = 0;
	density = 0.5;
	ecm_fiber_alignment.assign(3, 0.0);
}

ECM::ECM()
{
	ECM_Voxel template_ecm_voxel;

	// ECM_Cartesian_Mesh ecm_mesh;

	ecm_mesh.resize(1, 1, 1);

	make_ecm_units();

	// ecm_voxels.push_back(template_ecm_voxel);
}

void ECM::make_ecm_units(void)
{
	// you want to ... grab the centers and ordering from ecm_mesh.voxels ... one by one. This will make it First - resize the mesh, then Second, resize the ECM. And you can make that happen in the constructor also, which might be the best idea???

	for (int i = 0; i < ecm_mesh.voxels.size(); i++)
	{

		ecm_voxels.push_back(ecm_voxel);
		ecm_voxels[i].mesh_index = ecm_mesh.voxels[i].mesh_index;
	}

	initialize_ECM();
}

void ECM::resize_ecm_units_from_ecm_mesh(void)
{
	ecm_voxels.resize(0);

	for (int i = 0; i < ecm_mesh.voxels.size(); i++)
	{
		ecm_voxels.push_back(ecm_voxel);
		ecm_voxels[i].mesh_index = ecm_mesh.voxels[i].mesh_index;
		ecm_voxels[i].center = ecm_mesh.voxels[i].center;
		ecm_voxels[i].volume = ecm_mesh.voxels[i].volume;
	}

	initialize_ECM();
}

void ECM::initialize_ECM(void)
{

	if (ecm_mesh.voxels.size() != ecm_voxels.size())
	{
		std::cout << "Resize ECM mesh to match ECM voxels before initializing ECM units to initial values" << std::endl;
		std::cout << " hit Enter to continue:" << std::flush;
		std::cin.get();
	}
}

ECM_options::ECM_options()
{
}

void setup_extracellular_matrix(void)
{
	// DEPENDS ON MICROENVIRONMENT - CALL SETUP MICROENVIRONEMNT FIRST!!!!!

	ecm.ecm_mesh.resize(default_microenvironment_options.X_range[0], default_microenvironment_options.X_range[1],
						default_microenvironment_options.Y_range[0], default_microenvironment_options.Y_range[1], default_microenvironment_options.Z_range[0], default_microenvironment_options.Z_range[1],
						default_microenvironment_options.dx, default_microenvironment_options.dy, default_microenvironment_options.dz);
	ecm.resize_ecm_units_from_ecm_mesh();

	ecm.ecm_mesh.display_information(std::cout);

	if (PhysiCell::parameters.strings("ecm_orientation_setup") == "csv")
	{
		return initialize_ecm_from_csv();
	}

	// double initial_anisotropy = PhysiCell::parameters.doubles("initial_ecm_anisotropy");

	int density_ind = microenvironment.find_density_index("ecm_density");
	int anisotropy_ind = microenvironment.find_density_index("ecm_anisotropy");

	for (int n = 0; n < ecm.ecm_mesh.voxels.size(); n++)
	{
		ecm.ecm_voxels[n].density = microenvironment.density_vector(n)[density_ind];
		ecm.ecm_voxels[n].anisotropy = microenvironment.density_vector(n)[anisotropy_ind];

		// For random 2-D initalization
		if (PhysiCell::parameters.strings("ecm_orientation_setup") == "random")
		{
			double theta = 6.2831853071795864769252867665590 * uniform_random();
			ecm.ecm_voxels[n].ecm_fiber_alignment = {ecm.ecm_voxels[n].anisotropy * cos(theta), ecm.ecm_voxels[n].anisotropy * sin(theta), 0.0};
		}
		// for starburst initialization
		else if (PhysiCell::parameters.strings("ecm_orientation_setup") == "starburst")
		{
			std::vector<double> position = ecm.ecm_mesh.voxels[n].center;
			normalize(&position);
			ecm.ecm_voxels[n].ecm_fiber_alignment = {position[0], position[1], 0}; // oriented out (perpindeicular to concentric circles)
			normalize(&ecm.ecm_voxels[n].ecm_fiber_alignment);
		}
		// for circular initialization
		else if (PhysiCell::parameters.strings("ecm_orientation_setup") == "circular")
		{
			std::vector<double> position = ecm.ecm_mesh.voxels[n].center;
			normalize(&position);
			ecm.ecm_voxels[n].ecm_fiber_alignment = {position[1], -position[0], 0}; // oriented in cirlce
			normalize(&ecm.ecm_voxels[n].ecm_fiber_alignment);
		}

		else if (PhysiCell::parameters.strings("ecm_orientation_setup") == "horizontal")
		{
			ecm.ecm_voxels[n].ecm_fiber_alignment = {1.0, 0.0, 0.0};
			normalize(&ecm.ecm_voxels[n].ecm_fiber_alignment);
		}
		else if (PhysiCell::parameters.strings("ecm_orientation_setup") == "vertical")
		{
			ecm.ecm_voxels[n].ecm_fiber_alignment = {0.0, 1.0, 0.0};
			normalize(&ecm.ecm_voxels[n].ecm_fiber_alignment);
		}
		else if (PhysiCell::parameters.strings("ecm_orientation_setup") == "split")
		{
			std::vector<double> position = ecm.ecm_mesh.voxels[n].center;
			normalize(&position);

			if (position[1] <= 0)
			{
				ecm.ecm_voxels[n].ecm_fiber_alignment = {1, -1, 0}; // oriented out (perpindeicular to concentric circles)
				normalize(&ecm.ecm_voxels[n].ecm_fiber_alignment);
			}
			else
			{
				ecm.ecm_voxels[n].ecm_fiber_alignment = {1, 1, 0}; // oriented out (perpindeicular to concentric circles)
				normalize(&ecm.ecm_voxels[n].ecm_fiber_alignment);
			}
		}
		else
		{
			std::cout << "WARNING: NO ECM ORIENTATION SPECIFIED. FIX THIS!!!" << std::endl;
			std::cout << "Halting program!!!" << std::endl;
			abort();
			return;
		}

	}
}

void initialize_ecm_from_csv(void)
{
	pugi::xml_node node;

	node = xml_find_node(physicell_config_root, "microenvironment_setup");
	node = xml_find_node(node, "ecm_setup");

	bool ecm_setup_enabled = node.attribute("enabled").as_bool();
	if (!ecm_setup_enabled)
	{
		std::cout << "WARNING: ECM setup not enabled. But you are trying to initialize from a csv file. Fix this!!!" << std::endl;
		return;
	}
	std::string format = node.attribute("format").as_string();
	if (format != "csv")
	{
		std::cout << "ERROR: ECM setup format not recognized. Fix this!!!" << std::endl;
		std::cout << "Halting program!!!" << std::endl;
		exit(-1);
	}

	std::string csv_file = xml_get_string_value(node, "filename");
	std::string csv_folder = xml_get_string_value(node, "folder");

	// The .csv file needs to contain one row per voxel.
	// Each row is a vector of values as follows: [x coord, y coord, z coord, ecm_density, ecm_orientation_x, ecm_orientation_y]
	// Thus, your table should be of size #voxels x 6 (rows x columns)

	std::string filename = csv_folder + "/" + csv_file;

	// open file 
	std::ifstream file( filename, std::ios::in );
	if( !file )
	{ 
		std::cout << "ERROR: " << filename << " not found during ecm initialization. Quitting." << std::endl; 
		exit(-1);
	}

	// determine if header row exists 
	std::string line; 
	std::getline( file , line );
	char c = line.c_str()[0];
	if( c == 'X' || c == 'x' )
	{ 
		// do not support this with a header yet
		if ((line.c_str()[2] != 'Y' && line.c_str()[2] != 'y') || (line.c_str()[4] != 'Z' && line.c_str()[4] != 'z'))
		{
			std::cout << "ERROR: Header row starts with x but then not y,z? What is this? Exiting now." << std::endl;
			file.close();
			exit(-1);
		}
		std::vector< std::string> column_names; // this will include x,y,z (so make sure to skip those below)
		std::stringstream stream(line);
		std::string field;

		while (std::getline(stream, field, ','))
		{
			column_names.push_back(field);
		}

		std::vector<std::string> expected_column_names = {"x", "y", "z", "ecm_density", "ecm_orientation_x", "ecm_orientation_y"};
		if (column_names.size() != expected_column_names.size())
		{
			std::cout << "ERROR: Wrong number of columns in the header row of the .csv file specifying BioFVM initial conditions." << std::endl
					  << "\tExpected: " << expected_column_names.size() << std::endl
					  << "\tFound: " << column_names.size() << std::endl
					  << "\tRemember, your table should have dimensions #voxels x 6." << std::endl
					  << "\tThe header row (if present) should be: x,y,z,ecm_density,ecm_orientation_x,ecm_orientation_y" << std::endl;
			exit(-1);
		}
		for (int i = 0; i<column_names.size(); i++) // skip x,y,z by starting at 3, not 0
		{
			if (column_names[i] != expected_column_names[i])
			{
				std::cout << "ERROR: Wrong column name in the header row of the .csv file specifying BioFVM initial conditions." << std::endl
						  << "\tExpected: " << expected_column_names[i] << std::endl
						  << "\tFound: " << column_names[i] << std::endl
						  << "\tRemember, your table should have dimensions #voxels x 6." << std::endl
						  << "\tThe header row (if present) should be: x,y,z,ecm_density,ecm_orientation_x,ecm_orientation_y" << std::endl;
				exit(-1);
			}
		}
	}

	std::cout << "Loading ecm initial conditions from CSV file " << filename << " ... " << std::endl;
	std::vector<int> voxel_set = {}; // set to check that no voxel value is set twice

	while (std::getline(file, line))
	{
		std::vector<double> data;
		csv_to_vector(line.c_str(), data);

		if ((voxel_set.size() == 0) && (data.size() != 6))
		{
			std::cout << "WARNING: Wrong number of ecm values supplied in the .csv file specifying ecm initial conditions." << std::endl
					  << "\tExpected: 6" << std::endl
					  << "\tFound: " << data.size() << std::endl
					  << "\tRemember, save your csv with columns as: x, y, z, ecm_density, ecm_orientation_x, ecm_orientation_y" << std::endl;
		}

		std::vector<double> position = {data[0], data[1], data[2]};
		int voxel_ind = ecm.ecm_mesh.nearest_voxel_index(position);
		ecm.ecm_voxels[voxel_ind].density = data[3];
		ecm.ecm_voxels[voxel_ind].ecm_fiber_alignment = {data[4], data[5], 0.0};
		ecm.ecm_voxels[voxel_ind].anisotropy = norm(ecm.ecm_voxels[voxel_ind].ecm_fiber_alignment);
		
		for (unsigned int j = 0; j < voxel_set.size(); j++)
		{
			if (voxel_ind == voxel_set[j])
			{
				std::cout << "ERROR : the csv-supplied initial conditions for BioFVM repeat the same voxel. Fix the .csv file and try again." << std::endl
						  << "\tPosition that was repeated: " << position << std::endl;
				exit(-1);
			}
		}
		voxel_set.push_back(voxel_ind);
	}

	if (voxel_set.size() != ecm.ecm_voxels.size())
	{
		std::cout << "ERROR : Wrong number of voxels supplied in the .csv file specifying ecm initial conditions." << std::endl
				  << "\tExpected: " << ecm.ecm_voxels.size() << std::endl
				  << "\tFound: " << voxel_set.size() << std::endl
				  << "\tRemember, your table should have dimensions #voxels x 6." << std::endl;
		exit(-1);
	}

	file.close();

	return;
}

void copy_ecm_data_to_BioFVM(void)
{
	// This enables the use of rules to change the cell behaviors as well as aspects of visualization in the Studio
	// this is ONE WAY!!! DO NOT MODIFY WITHING BIOFVM!!! OR USE CELLS TO CHANGE THIS!!!!!!!!!

	// found it greatly increase run time - could possibly improve this by calling at the phenotype time step - thats when rules are applied
	// look into this later - for now just going the route of signals/behaviors or even raw Hill functions
	// Note also that if you want to use the ECM to update a mechanics thing (like speed), rules might not be the right approach anyway (called to infrequently) - it would depend on the concept/intend - is it a phenotypic change or a physical reaction to environment?

	static double next_copy_time = 0.0;
	if (PhysiCell_globals.current_time < next_copy_time - 0.5 * diffusion_dt)
	{
		return;
	}
	int number_of_voxels = ecm.ecm_mesh.voxels.size();

	static int ecm_anisotropy_index = microenvironment.find_density_index("ecm_anisotropy");
	if (ecm_anisotropy_index == -1)
	{
	    std::cout << "        static int ecm_anisotropy_index = " <<ecm_anisotropy_index << std::endl;
	    std::exit(-1);
	}
	static int ecm_density_index = microenvironment.find_density_index("ecm_density");
	if (ecm_density_index == -1)
	{
	    std::cout << "        static int ecm_density_index = " <<ecm_density_index << std::endl;
	    std::exit(-1);
	}
	static int ecm_orientation_x_index = microenvironment.find_density_index("ecm_orientation_x");
	if (ecm_orientation_x_index == -1)
	{
	    std::cout << "        static int ecm_orientation_x_index = " <<ecm_orientation_x_index << std::endl;
	    std::exit(-1);
	}
	static int ecm_orientation_y_index = microenvironment.find_density_index("ecm_orientation_y");
	if (ecm_orientation_y_index == -1)
	{
	    std::cout << "        static int ecm_orientation_y_index = " <<ecm_orientation_y_index << std::endl;
	    std::exit(-1);
	}
	// static int ecm_radial_index = microenvironment.find_density_index("ecm_radial");
	// if (ecm_radial_index == -1)
	// {
	//     std::cout << "        static int ecm_radial_index = " <<ecm_radial_index << std::endl;
	//     std::exit(-1);
	// }

	for (int n = 0; n < number_of_voxels; n++)
	{
		microenvironment.density_vector(n)[ecm_density_index] = ecm.ecm_voxels[n].density;
		microenvironment.density_vector(n)[ecm_anisotropy_index] = ecm.ecm_voxels[n].anisotropy;
		microenvironment.density_vector(n)[ecm_orientation_x_index] = ecm.ecm_voxels[n].ecm_fiber_alignment[0];
		microenvironment.density_vector(n)[ecm_orientation_y_index] = ecm.ecm_voxels[n].ecm_fiber_alignment[1];
		// if (ecm_density_index != -1)
		// {
		// 	microenvironment.density_vector(n)[ecm_density_index] = ecm.ecm_voxels[n].density;
		// }
		// if (ecm_anisotropy_index != -1)
		// {
		// 	microenvironment.density_vector(n)[ecm_anisotropy_index] = ecm.ecm_voxels[n].anisotropy;
		// }
		// if (ecm_orientation_x_index != -1)
		// {
		// 	microenvironment.density_vector(n)[ecm_orientation_x_index] = ecm.ecm_voxels[n].ecm_fiber_alignment[0];
		// }
		// if (ecm_orientation_y_index != -1)
		// {
		// 	microenvironment.density_vector(n)[ecm_orientation_y_index] = ecm.ecm_voxels[n].ecm_fiber_alignment[1];
		// }

		// std::vector<unsigned int> cart_inds = microenvironment.cartesian_indices(n);
		// std::vector<double> voxel_position = {microenvironment.mesh.x_coordinates[cart_inds[0]], microenvironment.mesh.y_coordinates[cart_inds[1]], microenvironment.mesh.z_coordinates[cart_inds[2]]};
		// if ((voxel_position[0] != 0 || voxel_position[1] != 0 || voxel_position[2] != 0) && ecm.ecm_voxels[n].anisotropy != 0)
		// {
		// 	microenvironment.density_vector(n)[ecm_radial_index] = fabs(dot_product(ecm.ecm_voxels[n].ecm_fiber_alignment, voxel_position)) / (ecm.ecm_voxels[n].anisotropy * norm(voxel_position));
		// }
		// else
		// {
		// 	microenvironment.density_vector(n)[ecm_radial_index] = 0;
		// }
	}

	next_copy_time = PhysiCell_globals.current_time + mechanics_dt; // mechanics_dt is the shortest time step at which the ecm is updated (movement==>update)
	return;
}
