////////
// title: PhyiCell/custom_modules/physicellmodule.cpp
// derived from: PhyiCell/main.cpp
//
// language: C/C++
// date: 2015-2024
// license: BSD-3-Clause
// author: Alexandre Bertin, Elmar Bucher, Paul Macklin
// original source code: https://github.com/MathCancer/PhysiCell
// modified source code: https://github.com/elmbeech/physicellembedding
// modified source code: https://github.com/Dante-Berth/PhysiGym
// input: https://docs.python.org/3/extending/extending.html
//
// description:
//   for the PhysiCell Python embedding the content of the regular main.cpp
//   was ported to this physicellmodule.cpp file.
////////


// load Python API
// since Python may define some pre-processor definitions which affect the standard headers on some systems, you must include Python.h before any standard headers are included.
#define PY_SSIZE_T_CLEAN
//#include <Python.h>

// load standard library
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <omp.h>

// loade PhysiCell library
#include "./core/PhysiCell.h"
#include "modules/PhysiCell_standard_modules.h"
#include "./custom_modules/custom.h"

// load namespace
using namespace BioFVM;
using namespace PhysiCell;

// global variable
char filename[1024];
std::ofstream report_file;
//std::vector<std::string> (*cell_coloring_function)(Cell*) = my_coloring_function;




int main( int argc, char* argv[] )
{
	// load and parse settings file(s)
	
	bool XML_status = false; 
	char copy_command [1024]; 
	if( argc > 1 )
	{
		XML_status = load_PhysiCell_config_file( argv[1] ); 
		sprintf( copy_command , "cp %s %s" , argv[1] , PhysiCell_settings.folder.c_str() ); 
	}
	else
	{
		XML_status = load_PhysiCell_config_file( "./config/PhysiCell_settings.xml" );
		sprintf( copy_command , "cp ./config/PhysiCell_settings.xml %s" , PhysiCell_settings.folder.c_str() ); 
	}
	if( !XML_status )
	{ exit(-1); }
	
	// copy config file to output directry 
	system( copy_command ); 
	
	// OpenMP setup
	omp_set_num_threads(PhysiCell_settings.omp_num_threads);
	
	// time setup 
	std::string time_units = "min"; 

     // reset global variables
    PhysiCell_globals = PhysiCell_Globals();  // bue 20240624: reset

    // Microenvironment setup //
    setup_microenvironment();  // modify this in the custom code

    // PhysiCell setup ///

    // set mechanics voxel size, and match the data structure to BioFVM
    double mechanics_voxel_size = 30;
    Cell_Container* cell_container = create_cell_container_for_microenvironment(microenvironment, mechanics_voxel_size);

    // Users typically start modifying here.
    random_seed();
    generate_cell_types();  // bue 20240624: delete cells; (re)load cell definitions;
    setup_tissue();
    // Users typically stop modifying here.

    // set MultiCellDS save options
    set_save_biofvm_mesh_as_matlab(true);
    set_save_biofvm_data_as_matlab(true);
    set_save_biofvm_cell_data(true);
    set_save_biofvm_cell_data_as_custom_matlab(true);

    // bue 20240624: reset mesh0
    BioFVM::reset_BioFVM_substrates_initialized_in_dom();

    // save initial data simulation snapshot
    //char filename[1024];  // bue 20240130: going global
    sprintf(filename, "%s/initial", PhysiCell_settings.folder.c_str());
    save_PhysiCell_to_MultiCellDS_v2(filename, microenvironment, PhysiCell_globals.current_time);

    // save data simulation snapshot output00000000
    if (PhysiCell_settings.enable_full_saves == true) {
        sprintf(filename, "%s/output%08u", PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index);
        save_PhysiCell_to_MultiCellDS_v2(filename, microenvironment, PhysiCell_globals.current_time);
    }

    // save initial svg cross section through z = 0 and legend
    PhysiCell_SVG_options.length_bar = 200;  // set cross section length bar to 200 microns
    std::vector<std::string> (*cell_coloring_function)(Cell*) = my_coloring_function;  // set pathology coloring function // bue 20240130: not going global

    sprintf(filename, "%s/legend.svg", PhysiCell_settings.folder.c_str());
    create_plot_legend(filename, cell_coloring_function);

    sprintf(filename, "%s/initial.svg", PhysiCell_settings.folder.c_str());
    SVG_plot(filename, microenvironment, 0.0, PhysiCell_globals.current_time, cell_coloring_function);

    // save svg cross section snapshot00000000
    if (PhysiCell_settings.enable_SVG_saves == true) {
        sprintf(filename, "%s/snapshot%08u.svg", PhysiCell_settings.folder.c_str(), PhysiCell_globals.SVG_output_index);
        SVG_plot(filename, microenvironment, 0.0, PhysiCell_globals.current_time, cell_coloring_function);
    }

    // save legacy simulation report
    //std::ofstream report_file;  // bue 20240130: going global
    if (PhysiCell_settings.enable_legacy_saves == true) {
        sprintf(filename, "%s/simulation_report.txt", PhysiCell_settings.folder.c_str());
        report_file.open(filename);  // create the data log file
        report_file << "simulated time\tnum cells\tnum division\tnum death\twall time" << std::endl;
        log_output(PhysiCell_globals.current_time, PhysiCell_globals.full_output_index, microenvironment, report_file);  // output00000000
    }

    // standard output
    display_citations();
    display_simulation_status(std::cout);  // output00000000

    // set the performance timers
    BioFVM::RUNTIME_TIC();
    BioFVM::TIC();
	
	// main loop 
	
	try 
	{		
		while( PhysiCell_globals.current_time < PhysiCell_settings.max_time + 0.1*diffusion_dt )
		{
			// save data if it's time. 
			if( fabs( PhysiCell_globals.current_time - PhysiCell_globals.next_full_save_time ) < 0.01 * diffusion_dt )
			{
				display_simulation_status( std::cout ); 
				if( PhysiCell_settings.enable_legacy_saves == true )
				{	
					log_output( PhysiCell_globals.current_time , PhysiCell_globals.full_output_index, microenvironment, report_file);
				}
				
				if( PhysiCell_settings.enable_full_saves == true )
				{	
					sprintf( filename , "%s/output%08u" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index ); 
					
					save_PhysiCell_to_MultiCellDS_v2( filename , microenvironment , PhysiCell_globals.current_time ); 
				}
				
				PhysiCell_globals.full_output_index++; 
				PhysiCell_globals.next_full_save_time += PhysiCell_settings.full_save_interval;
			}
			
			// save SVG plot if it's time
			if( fabs( PhysiCell_globals.current_time - PhysiCell_globals.next_SVG_save_time  ) < 0.01 * diffusion_dt )
			{
				if( PhysiCell_settings.enable_SVG_saves == true )
				{	
					sprintf( filename , "%s/snapshot%08u.svg" , PhysiCell_settings.folder.c_str() , PhysiCell_globals.SVG_output_index ); 
					SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function );
					
					PhysiCell_globals.SVG_output_index++; 
					PhysiCell_globals.next_SVG_save_time  += PhysiCell_settings.SVG_save_interval;
				}
			}

			// update the microenvironment
			microenvironment.simulate_diffusion_decay( diffusion_dt );
			
			// run PhysiCell 
			((Cell_Container *)microenvironment.agent_container)->update_all_cells( PhysiCell_globals.current_time );
			
			/*
			  Custom add-ons could potentially go here. 
			*/
			
			PhysiCell_globals.current_time += diffusion_dt;
		}
		
		if( PhysiCell_settings.enable_legacy_saves == true )
		{			
			log_output(PhysiCell_globals.current_time, PhysiCell_globals.full_output_index, microenvironment, report_file);
			report_file.close();
		}
	}
	catch( const std::exception& e )
	{ // reference to the base of a polymorphic object
		std::cout << e.what(); // information from length_error printed
	}
	
	// save a final simulation snapshot 
	
	sprintf( filename , "%s/final" , PhysiCell_settings.folder.c_str() ); 
	save_PhysiCell_to_MultiCellDS_v2( filename , microenvironment , PhysiCell_globals.current_time ); 
	
	sprintf( filename , "%s/final.svg" , PhysiCell_settings.folder.c_str() ); 
	SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function );
	
	// timer 
	
	std::cout << std::endl << "Total simulation runtime: " << std::endl; 
	BioFVM::display_stopwatch_value( std::cout , BioFVM::runtime_stopwatch_value() ); 


	return 0; 
}
