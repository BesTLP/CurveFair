#ifndef CURVEFAIR_H
#define CURVEFAIR_H

#include <Geom_BSplineCurve.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <Geom_BSplineCurve.hxx>
#include <GeomAPI_PointsToBSpline.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <GeomConvert_BSplineCurveToBezierCurve.hxx>
#include <GCPnts_AbscissaPoint.hxx>

#include <BRep_Tool.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <BRepPrim_FaceBuilder.hxx>

#include <TopoDS_Shape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>     
#include <TopExp_Explorer.hxx>

#include <vector>
#include <string>

#include <Eigen/Dense>
#include <math_GaussSingleIntegration.hxx>

#include <Interface_Static.hxx>
#include <math.hxx>
#include <mutex>
#include <map>
#include <math_Matrix.hxx>
#include "GeomConvert_CompCurveToBSplineCurve.hxx"
#include "SurfacemodelingTool.h"

const Standard_Real HAUSDORFFDISTANCETOL = 50;
const Standard_Real ALPHA = 1;
const Standard_Real ALPHARATIO = 0.8;
const Standard_Integer PARANUM = 100;
const Standard_Real FITTOLERANCE = 50;
class CurveFair
{
public:

    //! @brief 光顺算法
    //! @note 输入的 theBsplineCurve 和 theFitPoints 不能同时为空, 如果曲线为空, 则根据 theFitPoints 重新拟合曲线
    //! @param [in]theBSplineCurve 输入的曲线
    //! @param [in]theFitPoints 输入的拟合点
    //! @param [in]theHausdorffDistance 距离容差(默认值为50，若型值点非空，则光顺后曲线与型值点距离需在容差内; 否则，最终曲线与原始曲线的距高在容差内)
    CurveFair(
        Handle(Geom_BSplineCurve)& theBSplineCurve, 
        std::vector<gp_Pnt>& theFitPoints, 
        Standard_Real theHausdorffDistanceTol = HAUSDORFFDISTANCETOL)
    {
        if (theBSplineCurve.IsNull() && theFitPoints.size() == 0)
        {
            m_errorCode.push_back("CONSTRUCT ERROR::theBsplineCurve and theFitPoints cannot be both empty.");
            return;
        }

        this->m_OriginalCurve = theBSplineCurve;
        this->m_OriginalFitPoints = theFitPoints;
        this->m_ArcLengthMappingSampleNum = PARANUM;
        this->m_Alpha = ALPHA;
        this->m_AlphaRatio = ALPHARATIO;
        this->m_HausdorffDistanceTol = theHausdorffDistanceTol;
        if (m_OriginalCurve.IsNull() && m_OriginalFitPoints.size() > 0)
        {
            m_OriginalCurve = SampleAndFitBSpline(m_OriginalFitPoints, FITTOLERANCE);
        }

        // 插入更多节点向量, 增加自由度
        if (m_OriginalCurve->Knots().Size() <= 5)
        {
            m_OriginalCurve = InsertUniformKnotsInAllSpans(m_OriginalCurve, 3);
        }
        else if (m_OriginalCurve->Knots().Size() <= 10)
        {
            m_OriginalCurve = InsertUniformKnotsInAllSpans(m_OriginalCurve, 2);
        }
        else
        {
            m_OriginalCurve = InsertUniformKnotsInAllSpans(m_OriginalCurve, 1);
        }

        m_OriginalCurve = RefineCurveByCurvatureAuto(m_OriginalCurve, OccArrayConvertoVector(m_OriginalCurve->KnotSequence()), 1);
        this->m_FitPointParameters = ReCalculateFitPointParameters(theFitPoints, m_OriginalCurve);
        this->m_OriginalCurveKnotSequence = OccArrayConvertoVector(m_OriginalCurve->KnotSequence());
        this->m_OriginalCurveKnots = OccArrayConvertoVector(m_OriginalCurve->Knots());
        this->m_FirstPole = m_OriginalCurve->Pole(1);
        this->m_LastPole = m_OriginalCurve->Pole(m_OriginalCurve->NbPoles());
    
        Perform(m_OriginalCurve);
           
    }


    void UniformCurve(Handle(Geom_BSplineCurve)& curve);

    //! @brief 执行光顺算法
    //! @param [in]theCurve 输入的曲线
    void Perform(const Handle(Geom_BSplineCurve)& theCurve);

    //! @brief 迭代光顺算法
    //! @param [in]D0 原始控制点
    //! @param [in]D 新生成控制点
    void Iterator(Eigen::MatrixXd D0, Eigen::MatrixXd D);

    //! @brief 最近邻排序算法
    //! @param [in]points 输入的点集
    //! @return 排序后的点集
    static std::vector<gp_Pnt> ReorderPointsNearestNeighbor(const std::vector<gp_Pnt>& points);

    //! @brief 根据型值点拟合曲线
    //! @param [in]theOriginalCurve 输入的曲线
    //! @param [in]theSampleNum 采样点数
    //! @param [in]theFitPoints 采样的拟合点
    //! @param [in]theOccfitFlag 是否利用Occ型值点拟合，若为false，则使用自编算法
    //! @param [in]theMaxDegree 曲线最大度数
    //! @param [in]theContinuity 曲线连续性
    //! @param [in]theTolerance 拟合精度
    //! @return 拟合后的曲线
    static Handle(Geom_BSplineCurve) CurveFair::SampleAndFitBSpline(
        const Handle(Geom_BSplineCurve)& theOriginalCurve,
        Standard_Integer theSampleNum,
        std::vector<gp_Pnt>& theFitPoints,
        Standard_Boolean theOccfitFlag = Standard_True,
        Standard_Integer theMaxDegree = 3,
        GeomAbs_Shape theContinuity = GeomAbs_G2,
        Standard_Real theTolerance = FITTOLERANCE);

    //! @brief 根据型值点拟合曲线
    //! @param [in]theFitPoints 输入的型值点
    //! @param [in]theOccfitFlag 是否利用Occ型值点拟合，若为false，则使用自编算法
    //! @param [in]theResortFlag 是否重新排序型值点，若为true，则重新排序型值点
    //! @param [in]theMaxDegree 曲线最大度数
    //! @param [in]theContinuity 曲线连续性
    //! @param [in]theTolerance 拟合精度
    //! @return 拟合后的曲线
    static Handle(Geom_BSplineCurve) SampleAndFitBSpline(
        std::vector<gp_Pnt>& theFitPoints,
        Standard_Real theTolerance = FITTOLERANCE,
        Standard_Boolean theOccfitFlag = Standard_True,
        Standard_Boolean theResortFlag = Standard_True,
        Standard_Integer theMaxDegree = 3,
        GeomAbs_Shape theContinuity = GeomAbs_C3);


    Handle(Geom_BSplineCurve) RefineCurveByCurvatureAuto(
        const Handle(Geom_BSplineCurve)& theCurve,
        const std::vector<Standard_Real>& myKnotSeq,
        const Standard_Integer baseInsertNum);

    static Handle(Geom_BSplineCurve) InsertKnotsBetweenKnotSpan(
        const Handle(Geom_BSplineCurve)& theCurve,
        Standard_Real t1,
        Standard_Real t2,
        Standard_Integer nInsert);

    static Handle(Geom_BSplineCurve) InsertKnots(
        const Handle(Geom_BSplineCurve)& theCurve,
        const std::vector<Standard_Real>& params,
        Standard_Integer mult);

    //! @brief 在曲线中插入一个节点（参数值）
//! @param [in] theCurve  原始 B 样条曲线
//! @param [in] u         插入的参数值（节点）
//! @param [in] mult      插入的重数（默认 1）
//! @return 插入节点后的新曲线
    static Handle(Geom_BSplineCurve) InsertKnot(
        const Handle(Geom_BSplineCurve)& theCurve,
        Standard_Real u,
        Standard_Integer mult = 1);

    //! @brief 在每一对相邻节点区间内插入 nInsert 个均匀分布的新节点
    //! @param [in] theCurve     原始曲线
    //! @param [in] nInsert      每个区间插入的新节点数量
    //! @return 插值后的新曲线
    static Handle(Geom_BSplineCurve) InsertUniformKnotsInAllSpans(
        const Handle(Geom_BSplineCurve)& theCurve,
        Standard_Integer nInsert);

    //! @brief 获取弧长参数映射函数 f(s)
    //! @note 原始曲线C(t)不满足弧长参数化,计算弧长映射函数 使得C(f(s)) 参数 s 满足弧长参数化
    //! @param [in]theCurve 输入的曲线
    //! @param [in]theTolerance 拟合精度
    //! @return 弧长参数映射函数
    Handle(Geom_BSplineCurve) GetArclengthParameterMapping(
        const Handle(Geom_BSplineCurve)& theCurve, 
        const Standard_Real theTolerance = 1e-7);


    //! @brief 给定s，计算t = f(s)
    //! @param [in]theParameter 参数
    //! @param [in]k 阶数
    //! @return t = f(s)
    static Standard_Real f(const Standard_Real theParameter, const Standard_Integer k = 0);

    //! @brief 给定t，计算s = f逆(t)
    //! @param [in]theBSplineCurve 输入的曲线
    //! @param [in]t 参数
    //! @return t参数对应的弧长参数s
    Standard_Real f_inverse(const Handle(Geom_BSplineCurve)& theBSplineCurve, Standard_Real t);

    //! @brief 计算控制点权重矩阵
    //! @param [in]theCurve 输入的曲线
    //! @param [in]theFirPoints 输入的型值点
    //! @return 型值点对应在曲线上的弧长参数值
    static std::vector<std::pair<gp_Pnt, Standard_Real>> ReCalculateFitPointParameters(
        const std::vector<gp_Pnt>& theFitPoints,
        const Handle(Geom_BSplineCurve)& theCurve = nullptr);


    std::pair<std::vector<gp_Pnt>, std::vector <gp_Pnt>> SampleCurveWithArclengthMapping(
        const Handle(Geom_BSplineCurve)& theCurve,
        const Standard_Integer nSamples = 10);

    //! @brief 获取光顺曲线中间结果
    //! @param [in]theCurve 输入的曲线
    //! @param [in]M 能量矩阵
    //! @param [in]V 控制点权重矩阵
    //! @param [in]D0 原始控制点
    //! @param [out]D 新生成控制点
    //! @param [in]theAlpha 能量权重
    //! @return 光顺曲线中间结果
    Handle(Geom_BSplineCurve) GetTempFairCurve(const Handle(Geom_BSplineCurve)& theCurve,
        Eigen::MatrixXd M,
        Eigen::MatrixXd V,
        Eigen::MatrixXd& D0,
        Eigen::MatrixXd& D,
        Standard_Real theAlpha);

    //! @brief 计算光顺曲线带切向约束中间结果
    //! @param [in]theCurve 输入的曲线
    //! @param [in]M 能量矩阵
    //! @param [in]V 控制点权重矩阵
    //! @param [in]D0 原始控制点
    //! @param [out]D 新生成控制点
    //! @param [in]theAlpha 能量权重
    //! @return 光顺曲线中间结果
    Handle(Geom_BSplineCurve) CurveFair::GetTempFairCurveWithTangentConstraint(
        const Handle(Geom_BSplineCurve)& theCurve,
        Eigen::MatrixXd M,
        Eigen::MatrixXd V,
        Eigen::MatrixXd& D0,
        Eigen::MatrixXd& D,
        Standard_Real theAlpha);

    //! @note 利用连续性矩阵计算光顺曲线
    //! @param [in]theCurve 输入的曲线
    //! @param [in]M 能量矩阵
    //! @param [in]C 连续性矩阵
    //! @param [in]V 控制点权重矩阵
    //! @param [in]D0 原始控制点
    //! @param [out]D 新生成控制点
    //! @param [in]theAlpha 能量权重
    //! @return 光顺曲线中间结果
    Handle(Geom_BSplineCurve) GetTempContinuniousFairCurve(
        const Handle(Geom_BSplineCurve)& theCurve,
        Eigen::MatrixXd M,          // n×n 矩阵
        Eigen::MatrixXd C,          // m×n 矩阵（包含首尾控制点）
        Eigen::MatrixXd V,          // n×n 矩阵
        Eigen::MatrixXd& D0,        // n×3 初始控制点
        Eigen::MatrixXd& D,         // n×3 输出控制点
        Standard_Real alpha);

    //! @brief 计算能量矩阵
    //! @param [in]theBSplineCurve 输入的曲线
    //! @param [in]p 曲线阶数
    //! @param [in]tol 拟合精度
    //! @return 能量矩阵
    Eigen::MatrixXd ComputeEnergyMatrix(
        const Handle(Geom_BSplineCurve)& theBSplineCurve,
        const Standard_Integer p,
        const Standard_Real theTolerance = 1e-6);

    //! @brief 计算连续性矩阵
    //! @param [in]theBSplineCurve 输入的曲线
    //! @param [in]p 曲线阶数
    //! @return 连续性矩阵
    Eigen::MatrixXd ComputeContinuityMatrix(
        const Handle(Geom_BSplineCurve)& theBSplineCurve,
        const Standard_Integer p);


     //! @brief 根据端点切向约束计算约束矩阵 C 和右端项 h
     //! @param [in]theD0 原始控制点矩阵 (n x 3)
     //! @param [out]theH  输出参数，约束方程的右端项 h (4 x 1)
     //! @return Eigen::MatrixXd 约束矩阵 C (4 x 3*(n-2))
    Eigen::MatrixXd ComputeConstraintMatrix(
        const Eigen::MatrixXd& theD0,
        Eigen::VectorXd& theH);

    //! @brief 设置控制点权重矩阵
    //! @note 利用 GrevilleAbscissae 计算控制点权重
    //! @param [in]theCurve 输入的曲线
    //! @param [in]V 控制点权重矩阵
    void SetControlPointWeightMatrix(
        const Handle(Geom_BSplineCurve)& theCurve, 
        Eigen::MatrixXd& V);

    //! @brief 计算基函数导数
    //! @note  Nurbs Book 公式(2 - 9)
    //! @param [in]u 参数
    //! @param [in]i 控制点索引
    //! @param [in]p 曲线阶数
    //! @param [in]k 导数阶数
    //! @param [in]Knots 节点向量
    static Standard_Real BasisFunctionDerivative(
        const Standard_Real u,
        const Standard_Integer i,
        const Standard_Integer p,
        const Standard_Integer k,
        const std::vector<Standard_Real>& Knots);

    //! @brief 创建新的B样条曲线
    //! @param [in]theOriginalCurve 原始曲线
    //! @param [in]D 控制点
    //! @return 新的B样条曲线
    Handle(Geom_BSplineCurve) CreateNewBSplineCurve(
        const Handle(Geom_BSplineCurve)& theOriginalCurve, 
        const Eigen::MatrixXd& D);                 

    //! @brief 计算曲线与曲线之间的Hausdorff距离
    //! @param [in]theOriginalCurve 原始曲线
    //! @param [in]theOperateCurve 操作曲线
    //! @return 曲线与曲线之间的Hausdorff距离
    Standard_Real GetCurveCurveHausdorffDistance(
        const Handle(Geom_BSplineCurve) theOriginalCurve,
        const Handle(Geom_BSplineCurve) theOperateCurve);

    //! @brief 计算曲线与型值点之间的Hausdorff距离
    //! @param [in]theFitPointParams 型值点及其对应弧长参数
    //! @param [in]theOperateCurve 操作曲线
    //! @return 曲线与型值点之间的Hausdorff距离
    Standard_Real GetFitPointsCurveHausdorffDistance(
        const std::vector<std::pair<gp_Pnt, Standard_Real>> theFitPointParams,
        const Handle(Geom_BSplineCurve)& theOperateCurve);

    //! @brief 获取弧长参数映射函数
    //! @return 弧长参数映射函数
    inline Handle(Geom_BSplineCurve) GetArcLengthMappingFunction() { return m_ArcLengthMappingFunction; }

    //! @brief 获取光顺曲线
    //! @return 光顺曲线
    inline Handle(Geom_BSplineCurve) GetResult() { return m_ResultCurve; }

    //! @brief 将Occ数组转换为标准向量
    //! @param [in]theOccArray Occ数组
    inline static std::vector<Standard_Real> OccArrayConvertoVector(TColStd_Array1OfReal theOccArray)
    {
        std::vector<Standard_Real> VectorArray;
        for (Standard_Integer i = theOccArray.Lower(); i <= theOccArray.Upper(); i++)
        {
            VectorArray.push_back(theOccArray.Value(i));
        }
        return VectorArray;
    }

    //! @brief 将Occ数组转换为标准向量
    //! @param [in]theOccArray Occ数组
    inline static std::vector<gp_Pnt> OccArrayConvertoVector(TColgp_Array1OfPnt theOccArray)
    {
        std::vector<gp_Pnt> PntArray;
        for (int i = 1; i <= theOccArray.Size(); i++)
        {
            PntArray.push_back(theOccArray.Value(i));
        }
        return PntArray;
    }

    //! @brief 是否有错误信息
    //! @return 错误信息
    inline Standard_Boolean HasErrorMessage() const 
    {
        return m_errorCode.size() > 0;
    }

    //! @brief 获取错误信息
    //! @return 错误信息
    inline std::vector<std::string> GetErrorMessage() const 
    {
        return m_errorCode;
    }

    //! @brief 获取原始曲线
    //! @return 原始曲线
    inline Handle(Geom_BSplineCurve)& GetOriginalCurve() { return m_OriginalCurve; }

    //! @brief 获取原始型值点 
    //! @return 原始型值点
    inline std::vector<gp_Pnt>& GetFitPoints() { return m_OriginalFitPoints; }

    //! @brief 获取Hausdorff距离结果
    //! @return Hausdorff距离结果
    inline Standard_Real GetHausdorffDistanceResult() const { return m_HausdorffDistanceResult; }
private:

    // 弧长参数映射采样点数
    Standard_Real m_ArcLengthMappingSampleNum;

    // 弧长参数映射函数
    static Handle(Geom_BSplineCurve) m_ArcLengthMappingFunction;

    // 原始曲线
    Handle(Geom_BSplineCurve) m_OriginalCurve;
    // 原始型值点
    std::vector<gp_Pnt> m_OriginalFitPoints;
    // 原始控制点
    TColgp_Array1OfPnt m_OriginalPoles;
    // 原始节点向量
    std::vector<Standard_Real> m_OriginalCurveKnots;
    // 原始节点向量序列
    std::vector<Standard_Real> m_OriginalCurveKnotSequence;
    // 原始曲线首控制点
    gp_Pnt m_FirstPole;
    // 原始曲线尾控制点
    gp_Pnt m_LastPole;

    // 光顺曲线
    Handle(Geom_BSplineCurve) m_ResultCurve;

    // 如果需要分段光顺的话
    std::vector<Handle(Geom_BSplineCurve)> m_ResultCurveArray;

    // 能量矩阵
    Eigen::MatrixXd M;
    // 连续性矩阵
    Eigen::MatrixXd C;
    // 控制点权重矩阵
    Eigen::MatrixXd V;
    // 光顺能量权重
    Standard_Real m_Alpha;
    // 距离容差
    Standard_Real m_HausdorffDistanceTol;
    // 最终Hausdorff距离
    Standard_Real m_HausdorffDistanceResult;

    // 光顺迭代下降Ratio
    Standard_Real m_AlphaRatio;
    std::vector<std::pair<gp_Pnt, Standard_Real>> m_FitPointParameters;
private:
    std::vector<std::string> m_errorCode;

public:
    static std::string ExportFilePath;
};
#endif
