#ifndef CURVEFAIR_H
#define CURVEFAIR_H

#include<Geom_BSplineCurve.hxx>
#include<vector>
#include<string>
#include <Eigen/Dense>

const Standard_Real HAUSDORFFDISTANCETOL = 10;
const Standard_Real CONTROLPOINTOFFSETTOL = 100;
const Standard_Real ALPHA = 1e-3;
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

        TColStd_Array1OfReal KnotArray = theBSplineCurve->KnotSequence();
        for (Standard_Integer i = KnotArray.Lower(); i <= KnotArray.Upper(); i++)
        {
            myKnots.push_back(KnotArray.Value(i));
        }
        this->alpha = alpha;
        this->HausdorffDistanceTol = HausdorffDistanceTol;
        this->ControlPointOffsetTol = ControlPointOffsetTol; // 控制点偏差
        this->k = k; // 求 k 阶导数

        m_OriginalCurve = theBSplineCurve;
        FirstPole = theBSplineCurve->Pole(1); // 第一个控制点
        LastPole = theBSplineCurve->Pole(theBSplineCurve->NbPoles()); // 最后一个控制点

        Perform(m_OriginalCurve);

    }

    void Perform(const Handle(Geom_BSplineCurve)& theCurve);

    void Iterator(Eigen::MatrixXd D0, Eigen::MatrixXd D);

    Handle(Geom_BSplineCurve) GetCurrentFairCurve(const Handle(Geom_BSplineCurve)& theCurve, Eigen::MatrixXd M, Eigen::MatrixXd V, Eigen::MatrixXd D0, Eigen::MatrixXd& D, Standard_Real alpha);
    // 从文件中加载CurveFair
    static std::vector<CurveFair> GetCurveFairsFromFile(const std::string& filePath, const Standard_Real paraNum = 101);

    // 获取弧长参数化曲线
    Handle(Geom_BSplineCurve) ComputeInnerByArcReparam(const Handle(Geom_BSplineCurve)& theCurve, const Standard_Real tol = 1e-7);

    // 计算两点之间弧长
    Standard_Real GetLengthByParam(Standard_Real sParam, Standard_Real eParam, Standard_Real tol = 1e-7);

    // 计算三阶导
    static Standard_Real GetCurvatureDerivativeSquare(const Standard_Real, const Standard_Address);

    Standard_Real CurveFair::GetCurveCurvature(const Handle(Geom_BSplineCurve) theCurve, const Standard_Real t);
    // 计算F
    Standard_Real GetFByGaussIntegral(const Standard_Real tol = 1e-6);

    static Standard_Real OneBasicFun(
        const Standard_Real u,
        const Standard_Integer i,
        const Standard_Integer p,
        const std::vector<Standard_Real>& Knots);

    static Standard_Real BSplineBasis(
        Standard_Integer u,
        Standard_Integer i,
        Standard_Real p,
        const std::vector<Standard_Real>& knots);

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

    Eigen::MatrixXd ComputeEnergyMatrix(
        const Handle(Geom_BSplineCurve)& theBSplineCurve,
        const Standard_Integer p,
        const Standard_Integer k,
        const Standard_Real tol = 1e-6);

    Eigen::Matrix3d GetNewControlPoints(
        const Eigen::MatrixXd& M,       // 能量矩阵 M
        const Eigen::MatrixXd& V,       // 对角矩阵 V
        const Eigen::Matrix3d& D0,      // 原始控制点向量
        const Standard_Real alpha);

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
    static std::vector<Handle(Geom_BSplineCurve)> TempCurveArray;
    static std::string ExportFilePath;
    static Handle(Geom_BSplineCurve) inner;
private:
    Standard_Real paraNum;
    Handle(Geom_BSplineCurve) outer;
    std::vector<Standard_Real> myKnots;
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
};
#endif
