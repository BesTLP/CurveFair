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
#include <BRepTools.hxx>              // ���ڶ�ȡ BRep �ļ�
#include <BRep_Builder.hxx>           // ���ڹ��� BRep ����
#include <TopoDS_Shape.hxx>           // ���ڱ�ʾ������״
#include <TopoDS_Edge.hxx>            // ���ڱ�ʾ��
#include <TopExp_Explorer.hxx>        // ���ڱ���������״�ı�
#include <BRep_Tool.hxx>              // ������ȡ�ߵļ��α�ʾ
#include <TopoDS.hxx>
#include <GeomConvert_BSplineCurveToBezierCurve.hxx>
#include <Interface_Static.hxx>
const Standard_Real HAUSDORFFDISTANCETOL = 10;
const Standard_Real CONTROLPOINTOFFSETTOL = 100;
const Standard_Real ALPHA = 1;
class CurveFair
{
public:
    CurveFair(Handle(Geom_BSplineCurve) theBSplineCurve, Standard_Real paraNum, Standard_Integer k,
        Standard_Real alpha = ALPHA, Standard_Real ControlPointOffsetTol = CONTROLPOINTOFFSETTOL, Standard_Real HausdorffDistanceTol = HAUSDORFFDISTANCETOL,
        Standard_Boolean refit = Standard_False, Standard_Integer refitNums = 20)
    {
        outer = theBSplineCurve;
        std::vector<gp_Pnt> points;
        if (refit) // ���Ҫ������ϵĻ�
            outer = SampleAndFitBSpline(outer, refitNums, points);

        // Initialize
        this->paraNum = paraNum;
        this->fitParams = TColStd_Array1OfReal(outer->Knots().Lower(), outer->Knots().Lower() + paraNum - 1);
        this->fitPoints = TColgp_Array1OfPnt(outer->Knots().Lower(), outer->Knots().Lower() + paraNum - 1);
        this->myKnotSeq = ConvertToVector(outer->KnotSequence());
        this->myKnots = ConvertToVector(outer->Knots());
        this->alpha = alpha;
        this->HausdorffDistanceTol = HausdorffDistanceTol;
        this->ControlPointOffsetTol = ControlPointOffsetTol; // ���Ƶ�ƫ��
        this->k = k;
        this->m_OriginalCurve = outer;
        this->FirstPole = outer->Pole(1); // ��һ�����Ƶ�
        this->LastPole = outer->Pole(outer->NbPoles()); // ���һ�����Ƶ�

        CurveFair::TempInitialArray.clear();
        CurveFair::TempResultArray.clear();
        CurveFair::TempIteratorArray.clear();

        //ConvertBSplineToBezierSegments(m_OriginalCurve);
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

    // ����������µ�B��������
    static Handle(Geom_BSplineCurve) SampleAndFitBSpline(
        const Handle(Geom_BSplineCurve)& originalCurve,
        Standard_Integer numSamplePoints,
        std::vector<gp_Pnt>& points,
        Standard_Boolean useOccFit = Standard_False,
        Standard_Integer maxDegree = 3,
        GeomAbs_Shape continuity = GeomAbs_C3,
        Standard_Real tolerance = 1e-7);

    // ʹ�ò�ֵ�����������
    static Handle(Geom_BSplineCurve) SampleAndInterpolateBSpline(
        const Handle(Geom_BSplineCurve)& originalCurve,
        int numSamplePoints,
        bool periodic = false,
        GeomAbs_Shape continuity = GeomAbs_C3
    );
    std::vector<Handle(Geom_BezierCurve)> ConvertBSplineToBezierSegments(const Handle(Geom_BSplineCurve)& theCurve);

    // ���� BEZIER ����Ϊ STEP �ļ�
    void ExportBezierCurvesToSTEP(const std::vector<Handle(Geom_BezierCurve)>& bezierCurves, const std::string& filename);
    void SetControlPointWeightMatrix(const Handle(Geom_BSplineCurve)& theCurve, Eigen::MatrixXd& V);
    Standard_Real f_inverse(const Handle(Geom_BSplineCurve)& theBSplineCurve, Standard_Real t);

    void Perform(const Handle(Geom_BSplineCurve)& theCurve);

    void Iterator(Eigen::MatrixXd D0, Eigen::MatrixXd D);

    Handle(Geom_BSplineCurve) GetCurrentFairCurve(const Handle(Geom_BSplineCurve)& theCurve, 
        Eigen::MatrixXd M, 
        Eigen::MatrixXd C, 
        Eigen::MatrixXd V, 
        Eigen::MatrixXd& D0, 
        Eigen::MatrixXd& D, 
        Standard_Real alpha);

    // ��ȡ��������������
    Handle(Geom_BSplineCurve) ComputeInnerByArcReparam(const Handle(Geom_BSplineCurve)& theCurve, const Standard_Real tol = 1e-7);

    // ��������֮�仡��
    Standard_Real GetLengthByParam(Standard_Real sParam, Standard_Real eParam, Standard_Real tol = 1e-7);

    // �������׵�
    static Standard_Real GetCurvatureDerivativeSquare(const Standard_Real, const Standard_Address);

    // ����F
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

    Eigen::MatrixXd ComputeContinuityMatrix(
        const Handle(Geom_BSplineCurve)& theBSplineCurve,
        const Standard_Integer p,
        const Standard_Integer k);

    //  ��ȡ���µĿ��Ƶ��ԭʼ���Ƶ�֮��ľ���
    Standard_Real GetControlPointsOffset(
        const TColgp_Array1OfPnt theOriginalPoles,
        const TColgp_Array1OfPnt theOperatePoles);

    // �����µ� B ��������
    Handle(Geom_BSplineCurve) CreateNewBSplineCurve(
        const Handle(Geom_BSplineCurve)& originalCurve, // ԭʼ����
        const Eigen::MatrixXd& D);                 // �µĿ��Ƶ�

    //  ��ȡ�������ߺ�ԭʼ����֮���Hausdorff����(����)
    static Standard_Real GetCurveCurveHausdorffDistance(
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

    Standard_Real alpha; // ������ֵ
    Standard_Real HausdorffDistanceTol;
    Standard_Real ControlPointOffsetTol;
    Standard_Integer k; // k ����
    Eigen::MatrixXd M;
    Eigen::MatrixXd C;
    gp_Pnt FirstPole;
    gp_Pnt LastPole;
    Eigen::MatrixXd V;
    Standard_Real OriEnergy;
    Standard_Real FairEnergy;

};
#endif
