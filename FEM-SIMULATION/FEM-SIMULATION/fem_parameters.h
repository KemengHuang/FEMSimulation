#pragma once
#ifndef FEM_PARAMETERS_H
#define FEM_PARAMETERS_H

namespace FEM {
	const static double PI = 3.1415926535897932;
	const static double density = 0.01;
	const static double YoungModulus = 10.0;
	const static double PoissonRate = 0.49;
	const static double explicit_time_step = 0.0001;
	const static double implicit_time_step = 0.001;
	const static double lengthRate = YoungModulus / (2 * (1 + PoissonRate));
	const static double volumRate = YoungModulus* PoissonRate / ((1 + PoissonRate) * (1 - 2 * PoissonRate));
	const static double aniosScale = lengthRate;
}


#endif // !FEM_PARAMETERS_H

