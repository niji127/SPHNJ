// SPH solver 2023.6

#include "sim.cuh"
#include "function.cuh"
#include "cuda_macro.cuh"
#include "global.h"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <iostream>
#include <windows.h>

int main()
{
	try
	{
		// initialize
		Sim::parameter.initiate();
		Sim::parameter.showInfo();

		Sim *sim = new Sim;
		sim->initiate();
		sim->showMemoryUsed();

		Sim::data_pointer.initiate(sim->fluid_list,sim->solid_list,sim->virt_list,sim->cell_list);

		// gpu mode
		DeviceFunction::GPUCalulation(sim);

		// calculate
		Sleep(1000);
		delete sim;
		return 0;
	}
	catch (const std::exception &e)
	{
		std::cerr << e.what() << '\n';
		return -1;
	}

	std::cout << "no no no" << std::endl;
	return 0;
}