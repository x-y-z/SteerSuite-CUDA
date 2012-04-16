//
// Copyright (c) 2009-2010 Shawn Singh, Mubbasir Kapadia, Petros Faloutsos, Glenn Reinman
// See license.txt for complete license.
//

#ifndef __STEERLIB_ENERGY_EFFICIENCY_TECHNIQUE_H__
#define __STEERLIB_ENERGY_EFFICIENCY_TECHNIQUE_H__

/// @file EnergyEfficiencyTechnique.h
/// @brief Declares the SteerLib::EnergyEfficiencyBenchmarkTechnique class.

#include "Globals.h"
#include "interfaces/BenchmarkTechniqueInterface.h"
#include "util/GenericException.h"

namespace SteerLib {
	
	class STEERLIB_API EnergyEfficiencyBenchmarkTechnique : public SteerLib::BenchmarkTechniqueInterface
	{
		void init() {
			_benchmarkScoreComputed = false;
		}

		void update(SimulationMetricsCollector * simulationMetrics, float timeStamp, float dt) { 
			// nothing to do here
		}

		float getTotalBenchmarkScore(SimulationMetricsCollector * simulationMetrics) { 
			if (!_benchmarkScoreComputed) _computeTotalBenchmarkScore(simulationMetrics);
			return _totalBenchmarkScore;
		}

		float getAgentBenchmarkScore(unsigned int agentIndex, SimulationMetricsCollector * simulationMetrics) { 
			AgentMetrics * metrics = simulationMetrics->getAgentCollector(agentIndex)->getCurrentMetrics();
			return (metrics->integralOfKineticEnergy / metrics->totalTimeEnabled);
		}

		void printTotalScoreDetails(SimulationMetricsCollector * simulationMetrics, std::ostream & out) {
			if (!_benchmarkScoreComputed) _computeTotalBenchmarkScore(simulationMetrics);
			out << " average instantaneous energy over all agents: " << _totalBenchmarkScore << "\n";
		}

		void printAgentScoreDetails(unsigned int agentIndex, SimulationMetricsCollector * simulationMetrics, std::ostream & out) { 
			AgentMetrics * metrics = simulationMetrics->getAgentCollector(agentIndex)->getCurrentMetrics();
			out << " total (integral) energy over simulation: " << metrics->totalAcceleration << "\n";
			out << " total simulation time: " << metrics->totalTimeEnabled << "\n";
			out << " average instantaneous energy = total energy / simulation time = " << metrics->totalAcceleration / metrics->totalTimeEnabled << "\n";
		}

	protected:
		void _computeTotalBenchmarkScore(SimulationMetricsCollector * simulationMetrics) {
			float accumulator = 0.0f;
			for (unsigned int i=0; i < simulationMetrics->getNumAgents(); i++) {
				accumulator += getAgentBenchmarkScore(i, simulationMetrics);
			}
			_totalBenchmarkScore = accumulator / ((float)simulationMetrics->getNumAgents());
			_benchmarkScoreComputed = true;
		}

		float _totalBenchmarkScore;
		bool _benchmarkScoreComputed;
	};

} // end namespace SteerLib

#endif
