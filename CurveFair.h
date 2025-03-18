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
        this->ControlPointOffsetTol = ControlPointOffsetTol; // ���Ƶ�ƫ��
        this->k = k; // �� k �׵���

        m_OriginalCurve = theBSplineCurve;
        FirstPole = theBSplineCurve->Pole(1); // ��һ�����Ƶ�
        LastPole = theBSplineCurve->Pole(theBSplineCurve->NbPoles()); // ���һ�����Ƶ�

        Perform(m_OriginalCurve);

    }

    void Perform(const Handle(Geom_BSplineCurve)& theCurve);

    void Iterator(Eigen::MatrixXd D0, Eigen::MatrixXd D);

    Handle(Geom_BSplineCurve) GetCurrentFairCurve(const Handle(Geom_BSplineCurve)& theCurve, Eigen::MatrixXd M, Eigen::MatrixXd V, Eigen::MatrixXd D0, Eigen::MatrixXd& D, Standard_Real alpha);
    // ���ļ��м���CurveFair
    static std::vector<CurveFair> GetCurveFairsFromFile(const std::string& filePath, const Standard_Real paraNum = 101);

    // ��ȡ��������������
    Handle(Geom_BSplineCurve) ComputeInnerByArcReparam(const Handle(Geom_BSplineCurve)& theCurve, const Standard_Real tol = 1e-7);

    // ��������֮�仡��
    Standard_Real GetLengthByParam(Standard_Real sParam, Standard_Real eParam, Standard_Real tol = 1e-7);

    // �������׵�
    static Standard_Real GetCurvatureDerivativeSquare(const Standard_Real, const Standard_Address);

    Standard_Real CurveFair::GetCurveCurvature(const Handle(Geom_BSplineCurve) theCurve, const Standard_Real t);
    // ����F
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
        const Eigen::MatrixXd& M,       // �������� M
        const Eigen::MatrixXd& V,       // �ԽǾ��� V
        const Eigen::Matrix3d& D0,      // ԭʼ���Ƶ�����
        const Standard_Real alpha);

    //  ��ȡ���µĿ��Ƶ��ԭʼ���Ƶ�֮��ľ���
    Standard_Real GetControlPointsOffset(
        const TColgp_Array1OfPnt theOriginalPoles,
        const TColgp_Array1OfPnt theOperatePoles);

    // �����µ� B ��������
    Handle(Geom_BSplineCurve) CreateNewBSplineCurve(
        const Handle(Geom_BSplineCurve)& originalCurve, // ԭʼ����
        const Eigen::MatrixXd& D);                 // �µĿ��Ƶ�

    //  ��ȡ�������ߺ�ԭʼ����֮���Hausdorff����(����)
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

    Standard_Real alpha; // ������ֵ
    Standard_Real HausdorffDistanceTol;
    Standard_Real ControlPointOffsetTol;
    Standard_Integer k; // k ����
    Eigen::MatrixXd M;
    gp_Pnt FirstPole;
    gp_Pnt LastPole;
};
#endif
