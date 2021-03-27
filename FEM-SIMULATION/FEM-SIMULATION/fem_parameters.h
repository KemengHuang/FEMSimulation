#pragma once
#ifndef FEM_PARAMETERS_H
#define FEM_PARAMETERS_H

namespace FEM {
	double YoungModulus = 2;
	double PoissonRate = 0.4;
	double explicit_time_step = 0.001;
	double implicit_time_step = 0.01;
	double mass = 0.05;
	double lengthRate = YoungModulus / (2 * (1 + PoissonRate));
	double volumRate = YoungModulus * PoissonRate / ((1 + PoissonRate) * (1 - 2 * PoissonRate));
	double aniosScale = lengthRate;
}


#endif // !FEM_PARAMETERS_H

