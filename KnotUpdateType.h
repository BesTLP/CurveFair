#pragma once

//本文记录了在迭代拟合过程中，当前节点向量（控制点数目）无法满足规定容差时候如何在初始节点向量上更新当前节点
//
//MID_KNOT_BY_SINGLE_ERROR:遍历所有参数，找到误差最大的参数点。假设该参数点位于区间【ui，ui+1），则选取该区间中点为新节点
// 
//MID_KNOT_BY_INTERVAL_ERROR:遍历所有参数区间，找到误差最大的参数区间。选取该区间中点为新节点
// 
//PARAM_BASED_BY_INTERVAL_ERROR:遍历所有参数区间，找到误差最大的参数区间。已中间的参数点为基础计算新节点区间。
//reference:Deng, Chongyang, and Hongwei Lin. "Progressive and iterative approximation for least squares 
// B-spline curve and surface fitting." Computer-Aided Design 47 (2014): 32-44.
//

enum KONT_UPDATE_TYPE
{
	MID_KNOT_BY_SINGLE_ERROR,
	MID_KNOT_BY_INTERVAL_ERROR,
	PARAM_BASED_BY_INTERVAL_ERROR
};
