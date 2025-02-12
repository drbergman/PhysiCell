#include "librr_intracellular.h"

#include <sstream>
#include <iostream>

RoadRunnerIntracellular::RoadRunnerIntracellular() : Intracellular()
{
	intracellular_type = "sbml";
    std::cout << "====== " << __FUNCTION__ << "() intracellular_type=" << intracellular_type << std::endl;
    std::cout << "====== " << __FUNCTION__ << "() sbml_filename = " <<  sbml_filename << std::endl;
	parameters.clear();
}


// constructor using XML node
RoadRunnerIntracellular::RoadRunnerIntracellular(pugi::xml_node& node)
{
	intracellular_type = "roadrunner";
	initialize_intracellular_from_pugixml(node);
}

RoadRunnerIntracellular::RoadRunnerIntracellular(RoadRunnerIntracellular* copy) 
{
    update_time_step = copy->update_time_step;
    previous_update_time = PhysiCell::PhysiCell_globals.current_time;
    next_librr_run = PhysiCell::PhysiCell_globals.current_time + update_time_step;
	intracellular_type = copy->intracellular_type;
	sbml_filename = copy->sbml_filename;
	parameters = copy->parameters;
    delay_terms = copy->delay_terms;
}

void RoadRunnerIntracellular::initialize_intracellular_from_pugixml(pugi::xml_node& node)
{
    update_time_step = PhysiCell::intracellular_dt; // default to this, but overwrite below if defined in XML

	pugi::xml_node node_sbml = node.child( "sbml_filename" );
	if ( node_sbml )
	{ 
        sbml_filename = PhysiCell::xml_get_my_string_value (node_sbml); 
        std::cout << "\n------------- "  << __FUNCTION__ << ": sbml_filename = " << sbml_filename << std::endl;
    }

    pugi::xml_node node_update_time_step = node.child( "intracellular_dt" );
    if ( node_update_time_step )
    { 
        update_time_step = PhysiCell::xml_get_my_double_value (node_update_time_step); 
        std::cout << "\n------------- "  << __FUNCTION__ << ": intracellular_dt = " << update_time_step << std::endl;
    }
	
	pugi::xml_node node_species = node.child( "map" );
	while( node_species )
	{
        // ---------  substrate
        
		std::string substrate_name = node_species.attribute( "PC_substrate" ).value(); 
		if( substrate_name != "" )
		{
			std::string species_name = node_species.attribute( "sbml_species" ).value();
			substrate_species[substrate_name] = species_name;
            std::cout << "\n------------- "  << __FUNCTION__ << ": species_name= " << species_name << std::endl;
		}
        // ---------  custom_data
		std::string custom_data_name = node_species.attribute( "PC_custom_data" ).value(); 
		if( custom_data_name != "" )
		{
			std::string species_name = node_species.attribute( "sbml_species" ).value();
			custom_data_species[custom_data_name] = species_name;
		}
        
        
        // ---------  phenotype_data
        std::string phenotype_name = node_species.attribute( "PC_phenotype" ).value(); 
        
		if( phenotype_name != "" )
		{
			std::string species_name = node_species.attribute( "sbml_species" ).value();
			phenotype_species[phenotype_name] = species_name;
		}

		node_species = node_species.next_sibling( "map" ); 
	}
	
    std::cout << "  ------- substrate_species map:"  << std::endl;
    for(auto elm : substrate_species)
    {
        std::cout << "      "  << elm.first << " -> " << elm.second << std::endl;
    }
    std::cout << "  ------- custom_data_species map:"  << std::endl;
    for(auto elm : custom_data_species)
    {
        std::cout << "      "  << elm.first << " -> " << elm.second << std::endl;
    }
    std::cout << std::endl;

    std::cout << "  ------- phenotype_species map:"  << std::endl;
    for(auto elm : phenotype_species)
    {
        std::cout << "      "  << elm.first << " -> " << elm.second << std::endl;
    }
    std::cout << std::endl;

}


void RoadRunnerIntracellular::start()
{
    // called when a new cell is created; creates the unique 'rrHandle'
    rrc::RRVectorPtr vptr;

    rrHandle = createRRInstance();

    if ( !rrc::loadSBML(rrHandle, (sbml_filename).c_str() ) )
    {
        std::cerr << "------------->>>>>  Error while loading SBML file  <-------------\n\n";
        exit(-1);
    }

    int r = rrc::getNumberOfReactions(rrHandle);
    int m = rrc::getNumberOfFloatingSpecies(rrHandle);
    int b = rrc::getNumberOfBoundarySpecies(rrHandle);
    int p = rrc::getNumberOfGlobalParameters(rrHandle);
    int c = rrc::getNumberOfCompartments(rrHandle);

    std::string species_names_str = stringArrayToString(rrc::getFloatingSpeciesIds(rrHandle));
    std::stringstream iss(species_names_str);
    std::string species_name;
    int idx = 0;
    while (iss >> species_name)
    {
        species_result_column_index[species_name] = idx;
        idx++;
    }

    vptr = rrc::getFloatingSpeciesConcentrations(rrHandle);
    
    rrc::freeVector(vptr);
}

bool RoadRunnerIntracellular::need_update()
{
    return PhysiCell::PhysiCell_globals.current_time >= this->next_librr_run - 0.5 * PhysiCell::diffusion_dt;
}

// solve the intracellular model
void RoadRunnerIntracellular::update()
{
    static double start_time = 0.0;
    static int num_vals = 2; // start time and end time

    rrc::freeRRCData (this->result);

    this->result = rrc::simulateEx (this->rrHandle, start_time, PhysiCell::PhysiCell_globals.current_time - previous_update_time, num_vals);  // start time, end time, and number of points
    this->previous_update_time = PhysiCell::PhysiCell_globals.current_time;
    this->next_librr_run = PhysiCell::PhysiCell_globals.current_time + update_time_step;
}

double RoadRunnerIntracellular::get_parameter_value(std::string param_name)
{
    rrc::RRVectorPtr vptr;

    vptr = rrc::getFloatingSpeciesConcentrations(this->rrHandle);

    int offset = species_result_column_index[param_name];
    double res = vptr->Data[offset];
    rrc::freeVector(vptr);
    return res;
}
	
// rwh: might consider doing a multi-[species_name, value] "set" method
void RoadRunnerIntracellular::set_parameter_value(std::string species_name, double value)
{
    rrc::RRVectorPtr vptr;

    vptr = rrc::getFloatingSpeciesConcentrations(this->rrHandle);
    int idx = species_result_column_index[species_name];
    vptr->Data[idx] = value;
    rrc::setFloatingSpeciesConcentrations(this->rrHandle, vptr);
    rrc::freeVector(vptr);
}

RoadRunnerIntracellular* getRoadRunnerModel(PhysiCell::Phenotype& phenotype) {
	return static_cast<RoadRunnerIntracellular*>(phenotype.intracellular);
}

void RoadRunnerIntracellular::save_libRR(std::string path, std::string index)
{
	std::string state_file_name = path + "/states_" + index + ".dat";
	std::ofstream state_file( state_file_name );
	state_file << "---------  dummy output from save_libRR  ---------" << std::endl;
	state_file << "ID,state" << std::endl;
	for( auto cell : *PhysiCell::all_cells )
		state_file << cell->ID << "," << cell->phenotype.intracellular->get_state() << std::endl;
	state_file.close();
}

std::string RoadRunnerIntracellular::get_state()
{
    return sbml_filename;
}


int RoadRunnerIntracellular::update_phenotype_parameters(PhysiCell::Phenotype& phenotype)
{
    for(auto elm : phenotype_species)
    {
        // motility params
        if (elm.first[0] == 'm')
        {
            if (elm.first == "mms")
            {
                phenotype.motility.migration_speed = phenotype.intracellular->get_parameter_value(elm.second);
            }
            else if (elm.first == "mpt")
            {
                phenotype.motility.persistence_time = phenotype.intracellular->get_parameter_value(elm.second);
            }
            else if (elm.first == "mmb")
            {
                phenotype.motility.migration_bias = phenotype.intracellular->get_parameter_value(elm.second);
            }
            else
            {
            }
        }
        // death params
        else if (elm.first[0] == 'd')
        {
            if (elm.first == "da")
            {                
                phenotype.death.rates[0] = phenotype.intracellular->get_parameter_value(elm.second);
            }
            else if (elm.first == "dn")
            {
                phenotype.death.rates[1] = phenotype.intracellular->get_parameter_value(elm.second);
            }
            else
            {
            }
        }
        // secretion params
        else if (elm.first[0] == 's')
        {
            // parsing attribute and getting substrate name
            std::string s = elm.first;
            std::string delimiter = "_";

            size_t pos = 0;
            std::string token;
            while ((pos = s.find(delimiter)) != std::string::npos) {
                token = s.substr(0, pos);
                s.erase(0, pos + delimiter.length());
            }
            int sub_index = microenvironment.find_density_index(s);

            //transport types
            //uptake rate
            if (elm.first.substr(0,3) == "sur")
            {
                phenotype.secretion.uptake_rates[1] = phenotype.intracellular->get_parameter_value(elm.second);
            }
            //secretion rate
            else if (elm.first.substr(0,3) == "ssr")
            {
                phenotype.secretion.secretion_rates[sub_index] = phenotype.intracellular->get_parameter_value(elm.second);
            }
            //secretion density
            else if (elm.first.substr(0,3) == "ssd")
            {
                phenotype.secretion.saturation_densities[sub_index] = phenotype.intracellular->get_parameter_value(elm.second);
            }
            //net export rate
            else if (elm.first.substr(0,3) == "ser")
            {
                phenotype.secretion.net_export_rates[sub_index] = phenotype.intracellular->get_parameter_value(elm.second);
            }
            else
            {
            }
        }
        
        // cycle params
        else if (elm.first[0] == 'c')
        {
            if (elm.first.substr(0,3) == "ctr")
            {
                // parsing attribute and getting substrate name
                std::string s = elm.first;
                std::string delimiter = "_";

                size_t pos = 0;
                std::string token;
                int counter = 0;
                int start_index;
                while ((pos = s.find(delimiter)) != std::string::npos) {
                    token = s.substr(0, pos);
                    if (counter == 1)
                    {
                        start_index = atoi( token.c_str() );
                    }
                    s.erase(0, pos + delimiter.length());
                    counter += 1;
                }
                int end_index = atoi( s.c_str() );
                phenotype.cycle.data.transition_rate(start_index,end_index) = phenotype.intracellular->get_parameter_value(elm.second);
            }
            else
            {
            }
        }
        
        // volume params
        else if (elm.first[0] == 'v')
        {
            if (elm.first == "vtsc")
            {
                phenotype.volume.target_solid_cytoplasmic = phenotype.intracellular->get_parameter_value(elm.second);
            }
            else if (elm.first == "vtsn")
            {
                phenotype.volume.target_solid_nuclear = phenotype.intracellular->get_parameter_value(elm.second);
            }
            else if (elm.first == "vff")
            {
                phenotype.volume.target_fluid_fraction = phenotype.intracellular->get_parameter_value(elm.second);
            }
            else
            {
            }
        }
        else
        {
        }
        
    }
    return 0;
}


int RoadRunnerIntracellular::validate_PhysiCell_tokens(PhysiCell::Phenotype& phenotype)
{
    for(auto elm : phenotype_species)
    {
        // motility params
        if (elm.first[0] == 'm')
        {
            if (elm.first == "mms")
            {
            }
            else if (elm.first == "mpt")
            {
            }
            else if (elm.first == "mmb")
            {
            }
            else
            {
                std::cout<< std::endl;
                std::cout << "ERROR: There is no specified token parameters in the name of \"" << elm.first << "\" at motility parameters. Please take a look token specifications." << std::endl;
                std::cout<< std::endl;
                std::cout<< std::endl;
                exit(-1);
            }
        }
        // death params
        else if (elm.first[0] == 'd')
        {
            if (elm.first == "da")
            {                
            }
            else if (elm.first == "dn")
            {
            }
            else
            {
                std::cout<< std::endl;
                std::cout << "ERROR: There is no specified token parameters in the name of \"" << elm.first << "\" at death parameters. Please take a look token specifications." << std::endl;
                std::cout<< std::endl;
                std::cout<< std::endl;
                exit(-1);
            }
        }
        // secretion params
        else if (elm.first[0] == 's')
        {
            // parsing attribute and getting substrate name
            std::string s = elm.first;
            std::string delimiter = "_";
            size_t pos = 0;
            std::string token;
            while ((pos = s.find(delimiter)) != std::string::npos) {
                token = s.substr(0, pos);
                s.erase(0, pos + delimiter.length());
            }
            int sub_index = microenvironment.find_density_index(s);
            if ( sub_index < 0 )
            {
                std::cout<< std::endl;
                std::cout << "ERROR: There is no substrate named in the name of \"" << s << "\" at microenvironment. Please take a look token specifications." << std::endl;
                std::cout<< std::endl;
                std::cout<< std::endl;
                exit(-1);
            }
            
            if (elm.first.substr(0,3) == "sur")
            {
            }
            else if (elm.first.substr(0,3) == "ssr")
            {
            }
            else if (elm.first.substr(0,3) == "ssd")
            {
            }
            else if (elm.first.substr(0,3) == "ser")
            {
            }
            else
            {
                std::cout<< std::endl;
                std::cout << "ERROR: There is no specified token parameters in the name of \"" << elm.first << "\" at secretion parameters. Please take a look token specifications." << std::endl;
                std::cout<< std::endl;
                std::cout<< std::endl;
                exit(-1);
            }
        }
        else if (elm.first[0] == 'c')
        {
            if (elm.first.substr(0,3) == "ctr")
            {
                // getting num of phases
                int num_of_phases = (&(phenotype.cycle.model()))->phases.size();
                
                // getting start and end indices
                std::string s = elm.first;
                std::string delimiter = "_";
                size_t pos = 0;
                std::string token;
                int counter = 0;
                int start_index;
                while ((pos = s.find(delimiter)) != std::string::npos) {
                    token = s.substr(0, pos);
                    if (counter == 1)
                    {
                        start_index = atoi( token.c_str() );
                    }
                    s.erase(0, pos + delimiter.length());
                    counter += 1;
                }
                int end_index = atoi( s.c_str() );
                
                // validating the indices
                if ( start_index > num_of_phases - 1)
                {
                    std::cout<< std::endl;
                    std::cout << "ERROR: Given transition start index is beyond cycle indices. Please double check it." << std::endl;
                    std::cout<< std::endl;
                    std::cout<< std::endl;
                    exit(-1);
                }
                if ( end_index > num_of_phases - 1)
                {
                    std::cout<< std::endl;
                    std::cout << "ERROR: Given transition end index is beyond cycle indices. Please double check it." << std::endl;
                    std::cout<< std::endl;
                    std::cout<< std::endl;
                    exit(-1);
                }
            }
            else
            {
                std::cout<< std::endl;
                std::cout << "ERROR: There is no specified token parameters in the name of \"" << elm.first << "\" at cycle parameters. Please take a look token specifications." << std::endl;
                std::cout<< std::endl;
                std::cout<< std::endl;
                exit(-1);
            }
        }
        
        else if (elm.first[0] == 'v')
        {
            if (elm.first == "vtsc")
            {
            }
            else if (elm.first == "vtsn")
            {
            }
            else if (elm.first == "vff")
            {
            }
            else
            {
                std::cout<< std::endl;
                std::cout << "ERROR: There is no specified token parameters in the name of \"" << elm.first << "\" at volume parameters. Please take a look token specifications." << std::endl;
                std::cout<< std::endl;
                std::cout<< std::endl;
                exit(-1);
            }
        }
        else
        {
            std::cout<< std::endl;
            std::cout << "ERROR: There is no specified token parameters in the name of \"" << elm.first << "\" at phenotypic parameters. Please take a look token specifications." << std::endl;
            std::cout<< std::endl;
            std::cout<< std::endl;
            exit(-1);
        }
        
    }
    std::cout << "---- Specified PhysiCell tokens at config file are validated. ----- " << std::endl;
    
    return 0;
}

int RoadRunnerIntracellular::validate_SBML_species()
{
    // reading SBML
    rrHandle = createRRInstance();
    if ( !rrc::loadSBML(rrHandle, (sbml_filename).c_str() ) )
    {
        std::cerr << "------------->>>>>  Error while loading SBML file  <-------------\n\n";
        exit(-1);
    } 
    // getting Species Names
    std::string species_names_str = stringArrayToString(rrc::getFloatingSpeciesIds(rrHandle));
    std::stringstream iss(species_names_str);
    std::string species_name;
    
    std::vector<std::string> all_species {};
    
    int idx = 0;
    while (iss >> species_name)
    {
        species_result_column_index[species_name] = idx;
        all_species.push_back(species_name);
        idx++;
    }

    // Phenotype Species 
    for (auto elm : phenotype_species)
    {
        bool exist = 0;
        for (int i=0; i < all_species.size(); i++)
        {
            if ( all_species[i] == elm.second )
            {
               exist = 1; 
            }
            idx++;  
        }
        if (!exist)
        {
            std::cout<< std::endl;
            std::cout << "ERROR: The specified SBML species in the name of \"" << elm.second << "\" at phenotypic species. Please take a look SBML species specifications." << std::endl;
            std::cout<< std::endl;
            std::cout<< std::endl;
            exit(-1);
        }
    }
    
    // Substrate Species
    for (auto elm : substrate_species)
    {
        bool exist = 0;
        for (int i=0; i < all_species.size(); i++)
        {
            if ( all_species[i] == elm.second )
            {
               exist = 1; 
            }
            idx++;  
        }
        if (!exist)
        {
            std::cout<< std::endl;
            std::cout << "ERROR: The specified SBML species in the name of \"" << elm.second << "\" at substrate species. Please take a look SBML species specifications." << std::endl;
            std::cout<< std::endl;
            std::cout<< std::endl;
            exit(-1);
        }
    }    

    // custom data species
    for (auto elm : custom_data_species)
    {
        bool exist = 0;
        for (int i=0; i < all_species.size(); i++)
        {
            if ( all_species[i] == elm.second )
            {
               exist = 1; 
            }
            idx++;  
        }
        if (!exist)
        {
            std::cout<< std::endl;
            std::cout << "ERROR: The specified SBML species in the name of \"" << elm.second << "\" at substrate species. Please take a look SBML species specifications." << std::endl;
            std::cout<< std::endl;
            std::cout<< std::endl;
            exit(-1);
        }
    }    
    
    std::cout << "---- Specified SBML species at config file are validated. ----- " << std::endl;
    return 0;
}

int RoadRunnerIntracellular::create_custom_data_for_SBML(PhysiCell::Phenotype& phenotype)
{
    return 0; 
}
