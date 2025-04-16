#ifndef CURVEFAIR_H
#define CURVEFAIR_H

#include <Geom_BSplineCurve.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <vector>
#include <string>
#include <Eigen/Dense>
#include <GeomAPI_PointsToBSpline.hxx>
#include <Geom_BSplineCurve.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <math_GaussSingleIntegration.hxx>
#include <BRepPrim_FaceBuilder.hxx>
#include <BRepTools.hxx>              // 用于读取 BRep 文件
#include <BRep_Builder.hxx>           // 用于构建 BRep 数据
#include <TopoDS_Shape.hxx>           // 用于表示几何形状
#include <TopoDS_Edge.hxx>            // 用于表示边
#include <TopExp_Explorer.hxx>        // 用于遍历几何形状的边
#include <BRep_Tool.hxx>              // 用于提取边的几何表示
#include <TopoDS.hxx>

const Standard_Real HAUSDORFFDISTANCETOL = 10;
const Standard_Real CONTROLPOINTOFFSETTOL = 2;
const Standard_Real ALPHA = 1;
class CurveFair
{
public:
    CurveFair(Handle(Geom_BSplineCurve) theBSplineCurve, Standard_Real paraNum, Standard_Integer k,
        Standard_Real alpha = ALPHA, Standard_Real ControlPointOffsetTol = CONTROLPOINTOFFSETTOL, Standard_Real HausdorffDistanceTol = HAUSDORFFDISTANCETOL) :
        outer(theBSplineCurve),
        paraNum(paraNum),
        fitParams(outer->Knots().Lower(), outer->Knots().Lower() + paraNum - 1),
        fitPoints(outer->Knots().Lower(), outer->Knots().Lower() + paraNum - 1)
    {
        myKnotSeq = ConvertToVector(theBSplineCurve->KnotSequence());
        myKnots = ConvertToVector(theBSplineCurve->Knots());
        this->alpha = alpha;
        this->HausdorffDistanceTol = HausdorffDistanceTol;
        this->ControlPointOffsetTol = ControlPointOffsetTol; // 控制点偏差
        this->k = k;

        m_OriginalCurve = theBSplineCurve;
        FirstPole = theBSplineCurve->Pole(1); // 第一个控制点
        LastPole = theBSplineCurve->Pole(theBSplineCurve->NbPoles()); // 最后一个控制点

        CurveFair::TempInitialArray.clear();
        CurveFair::TempResultArray.clear();
        CurveFair::TempIteratorArray.clear();

        Perform(m_OriginalCurve);
    }

    std::vector<Standard_Real> ConvertToVector(TColStd_Array1OfReal OccArray)
    {
        std::vector<Standard_Real> VectorArray;
        for (Standard_Integer i = OccArray.Lower(); i <= OccArray.Upper(); i++)
        {
            VectorArray.push_back(OccArray.Value(i));
        }
        return VectorArray;
    }
    Standard_Real f_inverse(const Handle(Geom_BSplineCurve)& theBSplineCurve, Standard_Real t);

    void Perform(const Handle(Geom_BSplineCurve)& theCurve);

    void Iterator(Eigen::MatrixXd D0, Eigen::MatrixXd D);

    Handle(Geom_BSplineCurve) GetCurrentFairCurve(const Handle(Geom_BSplineCurve)& theCurve, Eigen::MatrixXd M, Eigen::MatrixXd V, Eigen::MatrixXd& D0, Eigen::MatrixXd& D, Standard_Real alpha);

    // 获取弧长参数化曲线
    Handle(Geom_BSplineCurve) ComputeInnerByArcReparam(const Handle(Geom_BSplineCurve)& theCurve, const Standard_Real tol = 1e-7);

    // 计算两点之间弧长
    Standard_Real GetLengthByParam(Standard_Real sParam, Standard_Real eParam, Standard_Real tol = 1e-7);

    // 计算三阶导
    static Standard_Real GetCurvatureDerivativeSquare(const Standard_Real, const Standard_Address);

    // 计算F
    Standard_Real GetFByGaussIntegral(const Standard_Real tol = 1e-6);

    static Standard_Real OneBasicFun(
        const Standard_Real u,
        const Standard_Integer i,
        const Standard_Integer p,
        const std::vector<Standard_Real>& Knots);

    static Standard_Real BasisFunctionDerivative(
        const Standard_Real u,
        const Standard_Integer i,
        const Standard_Integer p,
        const Standard_Integer k,
        const std::vector<Standard_Real>& Knots);

    static Standard_Real DerivativeSquaredTwoCallback(
        const Standard_Real u,
        const Standard_Address theAddress);

    static Standard_Real DerivativeSquaredThirdCallback(
        const Standard_Real u,
        const Standard_Address theAddress);

    Standard_Real ComputeDerivativeIntegral(
        const std::vector<Standard_Real> Knots,
        const Standard_Integer i,
        const Standard_Integer j,
        const Standard_Integer p,
        const Standard_Integer k,
        const Standard_Real tol = 1e-6);

    Standard_Real ComputeDerivativeIntegral(
        const std::vector<Standard_Real> Knots,
        const Standard_Integer i,
        const Standard_Integer p,
        const Standard_Integer k,
        const Standard_Real tol = 1e-6);

    Eigen::MatrixXd ComputeEnergyMatrix(
        const Handle(Geom_BSplineCurve)& theBSplineCurve,
        const Standard_Integer p,
        const Standard_Integer k,
        const Standard_Real tol = 1e-6);

    //  获取最新的控制点和原始控制点之间的距离
    Standard_Real GetControlPointsOffset(
        const TColgp_Array1OfPnt theOriginalPoles,
        const TColgp_Array1OfPnt theOperatePoles);

    // 创建新的 B 样条曲线
    Handle(Geom_BSplineCurve) CreateNewBSplineCurve(
        const Handle(Geom_BSplineCurve)& originalCurve, // 原始曲线
        const Eigen::MatrixXd& D);                 // 新的控制点

    //  获取最新曲线和原始曲线之间的Hausdorff距离(采样)
    Standard_Real GetCurveCurveHausdorffDistance(
        const Handle(Geom_BSplineCurve) theOriginalCurve,
        const Handle(Geom_BSplineCurve) theOperateCurve);

    static Standard_Real f(const Standard_Real theParameter, const Standard_Integer k = 0);

    inline Handle(Geom_BSplineCurve) GetOuter() { return outer; }
    inline Handle(Geom_BSplineCurve) GetInner() { return inner; }
    inline TColgp_Array1OfPnt GetFitPoint() { return fitPoints; }
    inline TColStd_Array1OfReal GetFitParam() { return fitParams; }

    inline Handle(Geom_BSplineCurve) GetResult() { return m_ResultCurve; }
public:
    static std::vector<Handle(Geom_BSplineCurve)> TempResultArray;
    static std::vector<Handle(Geom_BSplineCurve)> TempInitialArray;
    static std::vector<Handle(Geom_BSplineCurve)> TempIteratorArray;
    static std::string ExportFilePath;
    static Handle(Geom_BSplineCurve) inner;
private:
    Standard_Real paraNum;
    Handle(Geom_BSplineCurve) outer;
    std::vector<Standard_Real> myKnots;
    std::vector<Standard_Real> myKnotSeq;
    TColStd_Array1OfReal fitParams;
    TColgp_Array1OfPnt fitPoints;
    Handle(Geom_BSplineCurve) m_OriginalCurve;
    TColgp_Array1OfPnt m_OriginalPoles;
    Handle(Geom_BSplineCurve) m_ResultCurve;

    Standard_Real alpha; // 罚函数值
    Standard_Real HausdorffDistanceTol;
    Standard_Real ControlPointOffsetTol;
    Standard_Integer k; // k 次求导
    Eigen::MatrixXd M;
    gp_Pnt FirstPole;
    gp_Pnt LastPole;
    Eigen::MatrixXd V;
    Standard_Real OriEnergy;
    Standard_Real FairEnergy;

};
#endif
