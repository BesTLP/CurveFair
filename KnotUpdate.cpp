#include "KnotUpdate.h"
#include <algorithm>
#include <map>
#include <Interpolate.h>


KnotUpdate::KnotUpdate(Handle(Geom_BSplineCurve)& Bspline, const std::vector<double>& Sequences, const std::vector<gp_Pnt>& Pnts, const std::vector<double>& Params)
	:myParams(Params), myPnts(Pnts), bspline(Bspline), maxError(10000)
{
	std::for_each(Sequences.begin(), Sequences.end(), [&](double REALT) {myCurrentSequences.push_back(REALT); });
	updateKnotsAndMutis();
}

double KnotUpdate::SelfSingleUpdate(KONT_UPDATE_TYPE type)
{
	switch (type)
	{
	case MID_KNOT_BY_SINGLE_ERROR:
	{
		return selfUpdateForMidKnot();
		break;
	}
	case MID_KNOT_BY_INTERVAL_ERROR:
	{
		return selfUpdateForMidKnot_IntervalError();
		break;
	}
	case PARAM_BASED_BY_INTERVAL_ERROR:
	{
		return selfUpdateForLspia();
		break;
	}
	default:
		break;
	}
	return 0.0;
}

double KnotUpdate::selfUpdateForLspia()
{
	double maxSingleParamError, maxIntervalError, leftKnot, rightKnot, newKnot;
	maxSingleParamError = maxIntervalError = error(myParams[0], myPnts[0]);
	int maxSingleErrorParamIndex;
	int maxParamIntervalLeftIndex = 0;
	int maxParamIntervalRightIndex = 0;
	int knotIndex;
	int paramIndex = 0;
	//��ʼ�����ڵ�Ͳ����㣬Ҫ���������һ�������ڽڵ㡣ֻ���������һ���±��ǰһ��
	for (knotIndex = 0; knotIndex < myCurrentKnots.size() - 1; knotIndex++)
	{
		//��ȡ��ǰ�ڵ����� [leftKnot,rightKnot)
		leftKnot = myCurrentKnots[knotIndex];
		rightKnot = myCurrentKnots[knotIndex + 1];
		//��ȡλ�ڵ�ǰ�ڵ������ڵĲ����㣬����ͬ��������󵥵���������������
		double intervalError = 0; //��¼��ǰ�ڵ��������
		int count = 0;//��¼��ǰ������ٸ�����
		while (paramIndex < myParams.size())//ȷ����ǰ�����±���Ч
		{
			if (InterPolateTool::isGreaterThan(myParams[paramIndex], rightKnot))//��ǰ����������ǰ�����Ҷ˵�
			{
				break;
			}
			//��ǰ����λ��������
			double singleError = error(myParams[paramIndex], myPnts[paramIndex]);
			maxSingleParamError = (maxSingleParamError > singleError) ? maxSingleParamError : singleError;
			intervalError += singleError;
			paramIndex++;
			count++;
		}
		//�������±���Ч������������

		if (InterPolateTool::isGreaterThan(intervalError, maxIntervalError))//Ϊadmissible 䣬仯Ҫ��������
		{
			maxParamIntervalRightIndex = paramIndex - 1;
			maxParamIntervalLeftIndex = paramIndex - count;
			maxIntervalError = intervalError;
		}
	}

	//ȫ�������ϣ������½ڵ�
	newKnot = 0;
	for (int i = maxParamIntervalLeftIndex; i <= maxParamIntervalRightIndex; i++)
	{
		newKnot += myParams[i];
	}
	newKnot /= (maxParamIntervalRightIndex - maxParamIntervalLeftIndex + 1);
	//��newKnot���뵽��ǰ�ڵ㣬ע��Ƚ��ǲ��Ǻ����нڵ���ȣ���������Ҫ��������ظ����Ƿ񳬹�����
	auto index = checkNewKnot(newKnot);
	if (index != -1)
	{
		newKnot = myCurrentKnots[index] + (myCurrentKnots[index + 1] - myCurrentKnots[index]) / 20;
	}
	updateSequences(newKnot);
	updateKnotsAndMutis();
	maxError = maxSingleParamError;
	return newKnot;
}

double KnotUpdate::selfUpdateForMidKnot(bool isSingle)
{
	double maxSingleParamError, leftKnot, rightKnot, newKnot, maxParam;
	newKnot = 0.5;
	maxSingleParamError = error(myParams[0], myPnts[0]);
	int knotIndex, maxKnotIndex;
	int paramIndex = 0;
	//��ʼ�����ڵ�Ͳ����㣬Ҫ���������һ�������ڽڵ㡣ֻ���������һ���±��ǰһ��
	for (knotIndex = 0; knotIndex < myCurrentKnots.size() - 1; knotIndex++)
	{
		//��ȡ��ǰ�ڵ����� [leftKnot,rightKnot)
		leftKnot = myCurrentKnots[knotIndex];
		rightKnot = myCurrentKnots[knotIndex + 1];
		//��ȡλ�ڵ�ǰ�ڵ������ڵĲ����㣬����ͬ��������󵥵���������������
		int count = 0;//��¼��ǰ������ٸ�����
		while (paramIndex < myParams.size())//ȷ����ǰ�����±���Ч
		{
			if (InterPolateTool::isGreaterThan(myParams[paramIndex], rightKnot))//ǰǰҶ˵
			{
				break;
			}
			//��ǰ����λ��������
			double singleError = error(myParams[paramIndex], myPnts[paramIndex]);
			if (InterPolateTool::isGreaterThan(singleError, maxSingleParamError))
			{
				maxSingleParamError = singleError;
				maxKnotIndex = knotIndex;
				maxParam = myParams[paramIndex];
			}
			paramIndex++;
		}
	}

	if (isSingle)
	{
		if (InterPolateTool::isEqual(maxParam, 0))
		{
			newKnot = (myCurrentKnots[0] + myCurrentKnots[1]) / 2;
		}
		else if (InterPolateTool::isEqual(maxParam, 1))
		{
			newKnot = (myCurrentKnots[myCurrentKnots.size() - 1] + myCurrentKnots[myCurrentKnots.size() - 2]) / 2;
		}
		else
		{
			newKnot = (myCurrentKnots[maxKnotIndex] + myCurrentKnots[maxKnotIndex + 1]) / 2;
		}
	}
	updateSequences(newKnot);
	updateKnotsAndMutis();
	maxError = maxSingleParamError;
	return newKnot;
}

double KnotUpdate::selfUpdateForMidKnot_IntervalError()
{
	double maxSingleParamError, maxIntervalError, leftKnot, rightKnot, newKnot;
	maxSingleParamError = maxIntervalError = error(myParams[0], myPnts[0]);
	int knotIndex, maxKnotIndex;
	int paramIndex = 0;
	//��ʼ�����ڵ�Ͳ����㣬Ҫ���������һ�������ڽڵ㡣ֻ���������һ���±��ǰһ��
	for (knotIndex = 0; knotIndex < myCurrentKnots.size() - 1; knotIndex++)
	{
		//��ȡ��ǰ�ڵ����� [leftKnot,rightKnot)
		leftKnot = myCurrentKnots[knotIndex];
		rightKnot = myCurrentKnots[knotIndex + 1];
		//��ȡλ�ڵ�ǰ�ڵ������ڵĲ����㣬����ͬ��������󵥵���������������
		double intervalError = 0; //��¼��ǰ�ڵ��������
		while (paramIndex < myParams.size())//ȷ����ǰ�����±���Ч
		{
			if (InterPolateTool::isGreaterThan(myParams[paramIndex], rightKnot))//ǰǰҶ˵
			{
				break;
			}
			//ǰλ
			double singleError = error(myParams[paramIndex], myPnts[paramIndex]);
			maxSingleParamError = (maxSingleParamError > singleError) ? maxSingleParamError : singleError;
			intervalError += singleError;
			paramIndex++;
		}
		//�������±���Ч������������

		if (InterPolateTool::isGreaterThan(intervalError, maxIntervalError))//Ϊadmissible 䣬仯Ҫ
		{
			maxIntervalError = intervalError;
			maxKnotIndex = knotIndex;
		}
	}

	//ȫ�������ϣ������½ڵ�
	newKnot = 0;
	if (InterPolateTool::isEqual(maxKnotIndex, myCurrentKnots.size() - 1))
	{
		newKnot = (myCurrentKnots[myCurrentKnots.size() - 1] + myCurrentKnots[myCurrentKnots.size() - 2]) / 2;
	}
	else
	{
		newKnot = (myCurrentKnots[maxKnotIndex] + myCurrentKnots[maxKnotIndex + 1]) / 2;
	}
	updateSequences(newKnot);
	updateKnotsAndMutis();
	maxError = maxSingleParamError;
	return newKnot;
}

double KnotUpdate::error(double u, const gp_Pnt& P)
{
	return P.Distance(bspline->Value(u));
}

void KnotUpdate::updateKnotsAndMutis()
{
	if (myCurrentSequences.empty()) return;

	std::map<double, int> knotMap;

	// ʹ��map��ͳ��ÿ���ڵ���ظ�����
	for (double value : this->myCurrentSequences) {
		bool found = false;
		for (auto& knot : knotMap) {
			if (InterPolateTool::isEqual(value, knot.first)) {
				knot.second++;
				found = true;
				break;
			}
		}
		if (!found) {
			knotMap[value] = 1;
		}
	}
	myCurrentKnots.clear();
	myCurrentMutis.clear();
	// mapתƵknotsmultiplicities
	for (const auto& knot : knotMap) {
		myCurrentKnots.push_back(knot.first);
		myCurrentMutis.push_back(knot.second);
	}
}

void KnotUpdate::updateSequences()
{
	myCurrentSequences.clear();
	for (size_t i = 0; i < myCurrentKnots.size(); i++)
	{
		for (int j = 0; j < myCurrentMutis[i]; j++)
		{
			myCurrentSequences.push_back(myCurrentKnots[i]);
		}
	}
}

void KnotUpdate::updateSequences(double newKnot)
{
	myCurrentSequences.clear();
	bool notPush = true;
	for (size_t i = 0; i < myCurrentKnots.size(); i++)
	{

		for (int j = 0; j < myCurrentMutis[i]; j++)
		{
			myCurrentSequences.push_back(myCurrentKnots[i]);
		}
		if (notPush && InterPolateTool::isGreaterThan(newKnot, myCurrentKnots[i]) && InterPolateTool::isLessThan(newKnot, myCurrentKnots[i + 1]))
		{
			myCurrentSequences.push_back(newKnot);
			notPush = false;
		}
	}
}

int KnotUpdate::checkNewKnot(double knot)
{
	for (size_t i = 0; i < myCurrentKnots.size(); i++)
	{
		if (InterPolateTool::isEqual(knot, myCurrentKnots[i]))
		{
			return i;
		}
	}
	return -1;
}
