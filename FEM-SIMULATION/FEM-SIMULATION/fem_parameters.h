#pragma once
#ifndef FEM_PARAMETERS_H
#define FEM_PARAMETERS_H

namespace FEM {

const static double YoungModulus = 1;
const static double PoissonRate = 0.35;
const static double explicit_time_step = 0.001;
const static double implicit_time_step = 0.01;
const static double mass = 0.01;
const static double mass_inverse = 1.0 / mass;
const static double lengthRate = YoungModulus / (2 * (1 + PoissonRate));
const static double volumRate = YoungModulus * PoissonRate / ((1 + PoissonRate) * (1 - 2 * PoissonRate));
const static double aniosScale = lengthRate;

}


#endif // !FEM_PARAMETERS_H

