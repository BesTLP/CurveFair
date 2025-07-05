#pragma once
#include "RealCompare.h"
#include "KnotUpdateType.h"
#include <Geom_BSplineCurve.hxx>
#include <vector>
	//this class use to update knot in iterate approximate
	class KnotUpdate
	{
	public:
		KnotUpdate(Handle(Geom_BSplineCurve)& Bspline, const std::vector<double>& Sequences, const std::vector<gp_Pnt>& Pnts, const std::vector<double>& Params);

		double SelfSingleUpdate(KONT_UPDATE_TYPE type);

		double getMaxError()
		{
			return maxError;
		}

		std::vector<double> getSequences()
		{
			return myCurrentSequences;
		}

		~KnotUpdate() {}

	private:
		double selfUpdateForLspia();

		double selfUpdateForMidKnot(bool isSingle = true);

		double selfUpdateForMidKnot_IntervalError();
	private:

		double error(double u, const gp_Pnt& P);

		void updateKnotsAndMutis();

		void updateSequences();

		void updateSequences(double newKnot);

		int checkNewKnot(double knot);

	private:
		Handle(Geom_BSplineCurve)& bspline;

		std::vector<double> myCurrentSequences;
		std::vector<gp_Pnt> myPnts;
		std::vector<double> myParams;

		std::vector<double> myCurrentKnots;
		std::vector<int> myCurrentMutis;

		double maxError;
	};

