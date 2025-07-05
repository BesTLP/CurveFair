#pragma once

//���ļ�¼���ڵ�����Ϲ����У���ǰ�ڵ����������Ƶ���Ŀ���޷�����涨�ݲ�ʱ������ڳ�ʼ�ڵ������ϸ��µ�ǰ�ڵ�
//
//MID_KNOT_BY_SINGLE_ERROR:�������в������ҵ�������Ĳ����㡣����ò�����λ�����䡾ui��ui+1������ѡȡ�������е�Ϊ�½ڵ�
// 
//MID_KNOT_BY_INTERVAL_ERROR:�������в������䣬�ҵ�������Ĳ������䡣ѡȡ�������е�Ϊ�½ڵ�
// 
//PARAM_BASED_BY_INTERVAL_ERROR:�������в������䣬�ҵ�������Ĳ������䡣���м�Ĳ�����Ϊ���������½ڵ����䡣
//reference:Deng, Chongyang, and Hongwei Lin. "Progressive and iterative approximation for least squares 
// B-spline curve and surface fitting." Computer-Aided Design 47 (2014): 32-44.
//

enum KONT_UPDATE_TYPE
{
	MID_KNOT_BY_SINGLE_ERROR,
	MID_KNOT_BY_INTERVAL_ERROR,
	PARAM_BASED_BY_INTERVAL_ERROR
};
