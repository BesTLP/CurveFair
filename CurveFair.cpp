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

#include"curveFair.h"
#include "../../OpenCASCADE-7.7.0-vc14-64/opencascade-7.7.0/src/CPnts/CPnts_AbscissaPoint.cxx"
#include <RealCompare.h>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <SurfaceModelingTool.h>

// ���徲̬��Ա����
std::vector<Handle(Geom_BSplineCurve)> CurveFair::TempCurveArray;
Handle(Geom_BSplineCurve) CurveFair::inner;
std::string CurveFair::ExportFilePath;
std::string format_as(Eigen::MatrixXd M)
{
    std::stringstream ss;
    for (Standard_Integer i = 0; i < M.rows(); i++)
    {
        for (Standard_Integer j = 0; j < M.cols(); j++)
        {
            ss << M(i, j) << " ";
        }
        ss << std::endl;
    }
    return ss.str();
}
std::string format_as(TColgp_Array1OfPnt Poles)
{
    std::stringstream ss;
    ss << "[";
    for (Standard_Integer i = Poles.Lower(); i <= Poles.Upper(); i++)
    {
        ss << "(" << Poles.Value(i).X() << ", " << Poles.Value(i).Y() << ", " << Poles.Value(i).Z() << ")";
    }
    ss << "]" << std::endl;
    return ss.str();
}

std::string format_as(std::vector<Standard_Real> Knots)
{
    std::stringstream ss;
    ss << "[";
    for (Standard_Integer i = 0; i < Knots.size(); i++)
    {
        ss << Knots[i];
        if (i != Knots.size() - 1)
        {
            ss << ", ";
        }
    }
    ss << "]" << std::endl;
    return ss.str();
}
std::string format_as(const Handle(Geom_BSplineCurve) theCurve)
{
    std::stringstream ss;
    Standard_Real aFirstParameter = theCurve->FirstParameter();
    Standard_Real aLastParameter = theCurve->LastParameter();
    TColStd_Array1OfReal aKnotArray = theCurve->Knots();
    ss << "FirstParameter:" << aFirstParameter << ", " << std::endl;
    ss << "LastParameter:" << aLastParameter << std::endl;
    ss << "[";
    for (Standard_Integer i = 1; i <= aKnotArray.Size(); i++)
    {
        ss << aKnotArray.Value(i);
        if (i == aKnotArray.Size())
        {
            ss << "]";
        }
        else
        {
            ss << ",";
        }
    }
    return ss.str();
}

std::string format_as(const std::vector<Handle(Geom_BSplineCurve)> theCurveArray)
{
    std::stringstream ss;

    for (auto curve : theCurveArray)
    {
        ss << format_as(curve) << std::endl;
    }
    return ss.str();
}

//// ���ļ�������CurveFair����
//std::vector<CurveFair> CurveFair::GetCurveFairsFromFile(const std::string& filePath, const Standard_Real paraNum)
//{
//    // ��ȡ�ļ���׺
//    std::string extension = filePath.substr(filePath.find_last_of('.') + 1);
//    TopoDS_Shape boundary;
//    if (extension == "brep")
//    {
//        // ��ʼ���߽�Shape
//        BRep_Builder B1;
//        // ���ļ���ȡBRep����
//        BRepTools::Read(boundary, filePath.c_str(), B1);
//    }
//
//    // ����Shape�еı�
//    TopExp_Explorer explorer(boundary, TopAbs_EDGE);
//    std::vector<CurveFair> cfArray;
//    for (; explorer.More(); explorer.Next())
//    {
//        TopoDS_Edge edge = TopoDS::Edge(explorer.Current());
//
//        // ��ȡ�ߵļ��α�ʾ
//        TopLoc_Location loc;
//        Standard_Real first, last;
//        Handle(Geom_Curve) gcurve = BRep_Tool::Curve(edge, loc, first, last);
//        gcurve = Handle(Geom_Curve)::DownCast(gcurve->Copy());
//
//        Handle(Geom_BSplineCurve) aGeom_BSplineCurve = Handle(Geom_BSplineCurve)::DownCast(gcurve);
//        if (!aGeom_BSplineCurve.IsNull()) {
//            aGeom_BSplineCurve->Segment(first, last);
//            if (!aGeom_BSplineCurve.IsNull() && aGeom_BSplineCurve->IsKind(STANDARD_TYPE(Geom_BSplineCurve)))
//            {
//                CurveFair cf(aGeom_BSplineCurve, paraNum);
//                cfArray.push_back(cf);
//            }
//        }
//    }
//
//    return cfArray;
//}


/*
* �㷨: �� �����bspline���ܳ��ȣ��� ������ȳ���step =  bspline / knotNum  - 1���� ��ȡ���ȳ��ȶ�Ӧ��param��
*/
Standard_Real CurveFair::GetLengthByParam(Standard_Real sParam, Standard_Real eParam, Standard_Real tol)
{
    GeomAdaptor_Curve curveAdapter(outer);
    return GCPnts_AbscissaPoint::Length(curveAdapter, sParam, eParam, tol);
}

/*
* dervN: N�׵���1��2��3��
* bsplie����������
* ���أ�N�׵�
*/
std::string format_as(gp_Vec Vec)
{
    std::stringstream ss;
    ss << Vec.X() << ", " << Vec.Y() << ", " << Vec.Z() << ", Magnitude(): " << Vec.Magnitude();
    return ss.str();
}
Handle(Geom_BSplineCurve) CurveFair::ComputeInnerByArcReparam(const Handle(Geom_BSplineCurve)& theCurve, const Standard_Real tol)
{
    // bspline����
    Standard_Real firstParam = theCurve->FirstParameter();
    Standard_Real lastParam = theCurve->LastParameter();

    auto bsplineKnots = theCurve->Knots();
    const Standard_Integer LOWER = bsplineKnots.Lower();
    const Standard_Integer UPPER = bsplineKnots.Upper();

    // bspline: ���㻡��
    GeomAdaptor_Curve curveAdapter(theCurve);
    Standard_Real bsplineLen = GCPnts_AbscissaPoint::Length(curveAdapter, firstParam, lastParam, tol);
    Standard_Real avgLen = bsplineLen / (paraNum - 1.0);
    Standard_Real paramStep = (lastParam - firstParam) / (paraNum - 1.0);

    // ����������: �����µĽڵ�����
    TColStd_Array1OfReal reparamKnots(LOWER, LOWER + paraNum - 1);
    TColgp_Array1OfPnt reparamPoints(LOWER, LOWER + paraNum - 1);
    reparamKnots.SetValue(LOWER, firstParam);
    reparamPoints.SetValue(LOWER, gp_Pnt(firstParam, firstParam, firstParam));

    Standard_Real curLen = 0.0;
    Standard_Real curParam = firstParam;
    for (Standard_Integer i = LOWER + 1; i < paraNum; i++)
    {
        curLen += avgLen;

        curParam += paramStep;

        GCPnts_AbscissaPoint gap(curveAdapter, curLen, firstParam);
        Standard_Real tparam = gap.Parameter();

        reparamKnots.SetValue(i, curParam);
        reparamPoints.SetValue(i, gp_Pnt(tparam, tparam, tparam));
    }
    reparamKnots.SetValue(paraNum, lastParam);
    reparamPoints.SetValue(paraNum, gp_Pnt(lastParam, lastParam, lastParam));

    // Ϊÿһ���㣬��Ӳ���
    fitPoints = reparamPoints;
    /*fitParams = reparamKnots;*/
    inner = GeomAPI_PointsToBSpline(
        reparamPoints,
        reparamKnots,
        theCurve->Degree(),
        theCurve->Degree(),
        GeomAbs_C0,
        1e-7
    );

    return inner;
}


Standard_Real CurveFair::GetCurveCurvature(const Handle(Geom_BSplineCurve) theCurve, const Standard_Real t)
{
    return 0;
}
Standard_Real CurveFair::GetCurvatureDerivativeSquare(const Standard_Real t, const Standard_Address curF)
{
    CurveFair* cf = (CurveFair*)curF;
    Handle(Geom_BSplineCurve) inner = cf->GetInner();
    Handle(Geom_BSplineCurve) outer = cf->GetOuter();


    // inner�ĵ���
    gp_Pnt inPnt;
    gp_Vec inDerv1, inDerv2, inDerv3;
    inner->D3(t, inPnt, inDerv1, inDerv2, inDerv3);
    Standard_Real inParam = inPnt.X();
    Standard_Real inD1 = inDerv1.X();
    Standard_Real inD2 = inDerv2.X();
    Standard_Real inD3 = inDerv3.X();

    // outer�ĵ���
    gp_Pnt outPnt;
    gp_Vec outDerv1, outDerv2, outDerv3;
    outer->D3(inParam, outPnt, outDerv1, outDerv2, outDerv3);

    // D�ĵ���
    gp_Pnt d0 = outPnt;
    gp_Vec d1 = outDerv1 * inD1;
    //std::cout << "��λ������:" << d1.Magnitude() << std::endl;
    gp_Vec d2 = outDerv2 * inD1 * inD1 + outDerv1 * inD2;
    gp_Vec d3 = outDerv3 * pow(inD1, 3) + outDerv2 * (3 * inD1 * inD2) + outDerv1 * inD3;
    Standard_Real aCurvature = d2.Magnitude() / pow(d1.Magnitude(), 2);
    Standard_Real  aCurvatureDerivative = d3.Magnitude() / pow(d1.Magnitude(), 3);
    
    //std::cout << "�������������ʣ�������:" << aCurvature << std::endl;
    //std::cout << "�������������ʱ仯��(������)" << aCurvatureDerivative << std::endl;

    //// �������ʡ����ʱ仯��
    //Standard_Real aCurvature = d1.Crossed(d2).Magnitude() / pow(d1.Magnitude(), 3);
    //std::cout << "���淽������:" << aCurvature << std::endl;
    //gp_Vec dCurvatureDt = (d1.Crossed(d3) * d1.Magnitude() - d1.Crossed(d2) * (3.0 * d1.Dot(d2) / d1.Magnitude())) / pow(d1.Magnitude(), 4);
    //Standard_Real aCurvatureDerivative = dCurvatureDt.Magnitude() / d1.Magnitude();
    //std::cout << "���淽�����ʱ仯��" << aCurvatureDerivative << std::endl;
    return pow(aCurvatureDerivative, 2);
}



Standard_Real CurveFair::GetFByGaussIntegral(const Standard_Real tol)
{
    CPnts_MyGaussFunction gaussFunction;
    //POP pout WNT
    gaussFunction.Init(GetCurvatureDerivativeSquare, this);
      //FG.Init(f3d,(Standard_Address)&C);
    //Adaptor3d_Curve adaptor(this->GetOuter());
    // ʹ����������װ B-Spline ����
    Handle(GeomAdaptor_Curve) adaptor = new GeomAdaptor_Curve(this->GetOuter());

    std::cout << format_as(this->GetOuter()) << std::endl;
    std::cout << format_as(this->GetInner()) << std::endl;
    math_GaussSingleIntegration TheLength(
        gaussFunction, 
        this->GetOuter()->FirstParameter(),
        this->GetOuter()->LastParameter(),
        order(*adaptor), 
        tol
    );
    return Abs(TheLength.Value());
}

Eigen::MatrixXd CurveFair::ComputeEnergyMatrix(
    const Handle(Geom_BSplineCurve)& theBSplineCurve,
    const Standard_Integer p,
    const Standard_Integer k,
    const Standard_Real tol)
{
    Standard_Integer n = theBSplineCurve->NbPoles();
    TColgp_Array1OfPnt Poles = theBSplineCurve->Poles();
    std::vector<Standard_Real> Knots;
    for (Standard_Integer index = theBSplineCurve->KnotSequence().Lower(); index <= theBSplineCurve->KnotSequence().Upper(); index++)
        Knots.push_back(theBSplineCurve->KnotSequence().Value(index));

    Eigen::MatrixXd M(n, n);
    M.setZero();

    for (Standard_Integer i = 1; i < n - 1; ++i)
    {
        for (Standard_Integer j = i; j < n - 1; ++j)
        {
            Standard_Real integral = ComputeDerivativeIntegral(Knots, i, j, p, k, tol);
            M(i, j) = integral;

            // M(i, j)Ϊ�Գƾ���
            if (i != j)
            {
                M(j, i) = M(i, j);
            }
        }
    }
    return M;

    ////��֤�������󵼹�ʽ
    //Eigen::MatrixXd Basis(n, 1);
    //Basis.setZero();
    //Standard_Real FirstParameter = Knots[0];
    //Standard_Real LastParameter = Knots[Knots.size() - 1];
    //for (Standard_Real u = FirstParameter; u < LastParameter; u += (LastParameter - FirstParameter) / 50)
    //{
    //    for (Standard_Integer i = 0; i < n; i++)
    //    {
    //        Standard_Real BasisDerivate = BasisFunctionDerivative(u, i, p, 3, Knots);
    //        std::cout << "BasisDerivate:" << BasisDerivate << std::endl;
    //        Basis(i, 0) += BasisDerivate;
    //    }

    //    Eigen::MatrixXd Temp = D0.transpose() * Basis;
    //    std::cout << "D0" << std::endl << format_as(D0.transpose()) << std::endl;
    //    std::cout << "Basis" << std::endl << format_as(Basis) << std::endl;
    //    std::cout << "D0 * Basis" << std::endl << format_as(Temp) << std::endl;
    //    
    //    gp_Vec Derivate3 = gp_Vec(Temp(0, 0), Temp(1, 0), Temp(2, 0));
    //    std::cout << "Derivate3" << format_as(Derivate3) << std::endl;
    //}
}

struct BasisDerivativeIntegralContext
{
    Standard_Integer i;                   // ����������
    Standard_Integer p;                   // ����������
    Standard_Integer k;                   // ���������������׵���Ϊ3��
    std::vector<Standard_Real> Knots;     // �ڵ�����

    BasisDerivativeIntegralContext(
        Standard_Integer index_i,
        Standard_Integer degree,
        Standard_Integer derivOrder,
        const std::vector<Standard_Real>& knotVector)
        : i(index_i), p(degree), k(derivOrder), Knots(knotVector) {}
};


Standard_Real CurveFair::ComputeDerivativeIntegral(
    const std::vector<Standard_Real> Knots,
    const Standard_Integer i,
    const Standard_Integer j,
    const Standard_Integer p,
    const Standard_Integer k,
    const Standard_Real tol)
{
    // ������Ч����֤
    if (i < 0 || i + p + 1 >= Knots.size() || j < 0 || j + p + 1 >= Knots.size() || k > p) 
    {
        return 0.0;
    }

    // ȷ���������� - �������ľ�ȷ��������
    // ͨ��Ϊknots[i] -> kntos[i + p + 1]
    // ���������i��j�ķ��������ص�
    Standard_Real lower_i = Knots[i];
    Standard_Real upper_i = Knots[i + p + 1];
    Standard_Real lower_j = Knots[j];
    Standard_Real upper_j = Knots[j + p + 1];

    // ȷ����Ч��������
    Standard_Real lowerBound = std::max(lower_i, lower_j);
    Standard_Real upperBound = std::min(upper_i, upper_j);


    // ��������Ƿ���Ч
    if (std::abs(upperBound - lowerBound) < 1e-10)
    {
        return 0.0;
    }

    Standard_Integer Step = 200;
    Standard_Real StepParameter = (upperBound - lowerBound) / Step;
    Standard_Real AbsResult = 0;
    //for (Standard_Real parameter = lowerBound; parameter <= upperBound; parameter += StepParameter) 
    //{
    //    AbsResult += BasisFunctionDerivative(parameter, i, p, k, Knots) * BasisFunctionDerivative(parameter, j, p, k, Knots);
    //}
    //AbsResult *= StepParameter;


    for (Standard_Real parameter = lowerBound; parameter <= upperBound; parameter += StepParameter)
    {
        Standard_Real f0 = f(parameter);
        Standard_Real f1 = f(parameter, 1);
        Standard_Real f2 = f(parameter, 2);
        Standard_Real f3 = f(parameter, 3);

        Standard_Real First = BasisFunctionDerivative(f0, i, p, 3, Knots) * pow(f1, 3)
            + 3 * BasisFunctionDerivative(f0, i, p, 2, Knots) * f1 * f2
            + BasisFunctionDerivative(f0, i, p, 1, Knots) * f3;

        Standard_Real Second = BasisFunctionDerivative(f0, j, p, 3, Knots) * pow(f1, 3)
            + 3 * BasisFunctionDerivative(f0, j, p, 2, Knots) * f1 * f2
            + BasisFunctionDerivative(f0, j, p, 1, Knots) * f3;

        AbsResult += First * Second;
    }
    return AbsResult;

    //// ��������������
    //BasisDerivativeIntegralContext context(i, p, k, Knots);

    //// ������˹���ֺ����������ûص�
    //CPnts_MyGaussFunction gaussFunction;
    //gaussFunction.Init(DerivativeSquaredThirdCallback, &context);
    //// ִ�и�˹����
    //math_GaussSingleIntegration integrator(
    //    gaussFunction,
    //    lowerBound,
    //    upperBound,
    //    59,
    //    tol
    //);
    //return Abs(integrator.Value());
}
Standard_Real CurveFair::DerivativeSquaredTwoCallback(const Standard_Real u, const Standard_Address theAddress)
{
    // �ӵ�ַ�ָ������Ľṹ��
    BasisDerivativeIntegralContext* context =
        static_cast<BasisDerivativeIntegralContext*>(theAddress);

    Standard_Real fu = f(u);
    Standard_Real fu_1 = f(u, 1);
    Standard_Real fu_2 = f(u, 2);
    Standard_Integer i = context->i;
    Standard_Integer p = context->p;
    Standard_Integer DerivateK = context->k;
    std::vector<Standard_Real> Knots = context->Knots;


    Standard_Real Nj3_2 = BasisFunctionDerivative(fu, i, p, DerivateK, Knots);
    Standard_Real Nj3_1 = BasisFunctionDerivative(fu, i, p, DerivateK, Knots);

    Standard_Real first_1 = Nj3_2 * fu_1 * fu_1;
    Standard_Real first_2 = Nj3_1 * fu_2;
    Standard_Real Mi = first_1 + first_2;

    // ���ص���
    return Mi;
}

Standard_Real CurveFair::DerivativeSquaredThirdCallback(const Standard_Real u, const Standard_Address theAddress)
{
    // �ӵ�ַ�ָ������Ľṹ��
    BasisDerivativeIntegralContext* context =
        static_cast<BasisDerivativeIntegralContext*>(theAddress);

    Standard_Integer i = context->i;
    Standard_Integer p = context->p;
    Standard_Integer DerivateK = context->k;
    std::vector<Standard_Real> Knots = context->Knots;

    Standard_Real fu = f(u);
    Standard_Real fu_1 = f(u, 1);
    Standard_Real fu_2 = f(u, 2);
    Standard_Real fu_3 = f(u, 3);
    Standard_Real Nj3_3 = BasisFunctionDerivative(fu, i, p, DerivateK, Knots);
    Standard_Real Nj3_2 = BasisFunctionDerivative(fu, i, p, DerivateK, Knots);
    Standard_Real Nj3_1 = BasisFunctionDerivative(fu, i, p, DerivateK, Knots);

    Standard_Real Mi = Nj3_3 * pow(fu_1, 3) + 3 * Nj3_2 * fu_1 * fu_1 + Nj3_1 * fu_3;
    // ���ص���
    return Mi;
}

// To compute the value of a b-spline basic function value 
// write by: xin.c
Standard_Real CurveFair::OneBasicFun(const Standard_Real u, const Standard_Integer i, const Standard_Integer p, const std::vector<Standard_Real>& Knots)
{
    Standard_Real Nip, uleft, uright, saved, temp;
    Standard_Integer m = Knots.size() - 1;
    std::vector<Standard_Real>N(p + 1);
    if ((i == 0 && isEqual(u, Knots[0])) || (i == m - p - 1 && isEqual(u, Knots[m])))
    {
        return 1.0;
    }

    if (isLessThan(u, Knots[i]) || isGreaterThanOrEqual(u, Knots[i + p + 1]))
    {
        return 0.0;
    }


    for (size_t j = 0; j <= p; j++)
    {
        if (isGreaterThanOrEqual(u, Knots[i + j]) && isLessThan(u, Knots[i + j + 1]))
        {
            N[j] = 1.0;
        }
        else
        {
            N[j] = 0.0;
        }
    }
    for (size_t k = 1; k <= p; k++)
    {
        if (N[0] == 0.0)
        {
            saved = 0.0;
        }
        else
        {
            saved = ((u - Knots[i]) * N[0]) / (Knots[i + k] - Knots[i]);
        }
        for (size_t j = 0; j < p - k + 1; j++)
        {
            uleft = Knots[i + j + 1];
            uright = Knots[i + j + k + 1];
            if (N[j + 1] == 0.0)
            {
                N[j] = saved;
                saved = 0.0;
            }
            else
            {
                temp = N[j + 1] / (uright - uleft);
                N[j] = saved + (uright - u) * temp;
                saved = (u - uleft) * temp;
            }
        }
    }
    Nip = N[0];
    return Nip;
}
/**
 * ����B������������k�׵���
 * ���ڵݹ鹫ʽ: N^{(k)}_{i,p} = p * ((N^{(k-1)}_{i,p-1} / (u_{i+p} - u_i)) - (N^{(k-1)}_{i+1,p-1} / (u_{i+p+1} - u_{i+1})))
 * ���Բο� Nurbs Book ��ʽ(2-9)
 *
 * @param u ����ֵ
 * @param i ����������
 * @param p ����������
 * @param k �󵼽���
 * @param Knots �ڵ�����
 * @return k�׵���ֵ
 */

Standard_Real CurveFair::f(const Standard_Real theParameter, const Standard_Integer k)
{
    if (k == 0) return inner->Value(theParameter).X();
    gp_Vec DerivateVector = inner->DN(theParameter, k);
    return DerivateVector.X();
}

Standard_Real CurveFair::BasisFunctionDerivative(const Standard_Real u, const Standard_Integer i, const Standard_Integer p,
    const Standard_Integer k, const std::vector<Standard_Real>& Knots)
{
    // 0�׵���ֱ�ӷ��ػ�����ֵ
    if (k == 0)
    {
        return OneBasicFun(u, i, p, Knots);
    }

    // ��k > pʱ������Ϊ0, ��p = 0ʱ������Ϊ0
    if (k > p || p == 0)
    {
        return 0.0;
    }

    // �����һ��: (p / (u_{i+p} - u_i)) * N^{(k-1)}_{i,p-1}
    Standard_Real Term1 = 0.0;
    Standard_Real Denominator1 = Knots[i + p] - Knots[i];
    if (std::abs(Denominator1) > 1e-10)
    {
        Standard_Real coef1 = p / Denominator1;
        Standard_Real BasisDerivative = BasisFunctionDerivative(u, i, p - 1, k - 1, Knots);
        Term1 = coef1 * BasisDerivative;
    }

    // ����ڶ���: (p / (u_{i+p+1} - u_{i+1})) * N^{(k-1)}_{i+1,p-1}
    Standard_Real Term2 = 0.0;
    Standard_Real Denominator2 = Knots[i + p + 1] - Knots[i + 1];
    if (std::abs(Denominator2) > 1e-10)
    {
        Standard_Real coef2 = p / Denominator2;
        Standard_Real BasisDerivative = BasisFunctionDerivative(u, i + 1, p - 1, k - 1, Knots);
        Term2 = coef2 * BasisDerivative;
    }

    // ���ݹ�ʽ����k�׵���
    return Term1 - Term2;
}

Eigen::Matrix3d CurveFair::GetNewControlPoints(const Eigen::MatrixXd& M, const Eigen::MatrixXd& V, const Eigen::Matrix3d& D0, const Standard_Real alpha)
{
    return (alpha * M + V).inverse() * V * D0;
}

Handle(Geom_BSplineCurve) CurveFair::CreateNewBSplineCurve(const Handle(Geom_BSplineCurve)& originalCurve, const Eigen::MatrixXd& D)
{

    // ��ȡԭʼ���ߵ� Degree
    Standard_Integer degree = originalCurve->Degree();
    // ��ȡԭʼ���ߵ� �ڵ������ʹ���
    TColStd_Array1OfReal knots(1, originalCurve->NbKnots());
    originalCurve->Knots(knots);

    TColStd_Array1OfInteger multiplicities(1, originalCurve->NbKnots());
    originalCurve->Multiplicities(multiplicities);

    // �� Eigen::Matrix3d ת��Ϊ OCC ���Ƶ�����
    TColgp_Array1OfPnt Poles(1, D.rows()); // OCC ���Ƶ�����
    for (Standard_Integer i = 0; i < D.rows(); i++)
    {
        Poles.SetValue(i + 1, gp_Pnt(D(i, 0), D(i, 1), D(i, 2)));
    }
    // �����µ� B ��������
    Handle(Geom_BSplineCurve) newCurve = new Geom_BSplineCurve(
        Poles, knots, multiplicities, degree);

    return newCurve;
}

Standard_Real CurveFair::GetControlPointsOffset(const TColgp_Array1OfPnt theOriginalPoles, const TColgp_Array1OfPnt theOperatePoles)
{
    if (theOriginalPoles.Size() != theOperatePoles.Size())
    {
        return -1;
    }

    Standard_Real Offset = INT_MIN;
    for (Standard_Integer i = theOriginalPoles.Lower(); i <= theOriginalPoles.Upper(); i++)
    {
        Offset = std::max(Offset, theOriginalPoles.Value(i).Distance(theOperatePoles.Value(i)));
    }
    return Offset;
}
Standard_Real CurveFair::GetCurveCurveHausdorffDistance(const Handle(Geom_BSplineCurve) theOriginalCurve, const Handle(Geom_BSplineCurve) theOperateCurve)
{
    Standard_Real FirstParameter = theOriginalCurve->FirstParameter();
    Standard_Real LastParameter = theOriginalCurve->LastParameter();
    const Standard_Integer numPoints = 100;
    const Standard_Real step = (LastParameter - FirstParameter) / numPoints;
    Standard_Real HausdorffDistance = INT_MIN;
    for (Standard_Real t = FirstParameter; t <= LastParameter; t += step)
    {
        gp_Pnt originPoint = theOriginalCurve->Value(t);
        GeomAPI_ProjectPointOnCurve projector(originPoint, theOperateCurve);

        // ����Ƿ��ҵ�ͶӰ��
        if (!projector.NbPoints())
        {
            std::cerr << "No projection found for the point on the curve." << std::endl;
        }

        Standard_Real minDistance = INT_MAX;
        for (Standard_Integer i = 1; i <= projector.NbPoints(); i++)
        {
            minDistance = std::min(minDistance, projector.Distance(i));
        }
          
        HausdorffDistance = std::max(HausdorffDistance, minDistance);
    }
    return HausdorffDistance;
}

void CurveFair::Perform(const Handle(Geom_BSplineCurve)& theCurve)
{
    ComputeInnerByArcReparam(theCurve);
    M = ComputeEnergyMatrix(m_OriginalCurve, m_OriginalCurve->Degree(), k);
    TColgp_Array1OfPnt Poles = theCurve->Poles();
    Standard_Integer n = theCurve->NbPoles();
    Eigen::MatrixXd V(n, n);
    V.setIdentity();  // ��ʼ��VΪ��λ����

    Eigen::MatrixXd D0(n, 3);
    Eigen::MatrixXd D(n, 3);
    for (Standard_Integer i = 0; i < n; i++)
    {
        D0(i, 0) = Poles.Value(i + 1).X();
        D0(i, 1) = Poles.Value(i + 1).Y();
        D0(i, 2) = Poles.Value(i + 1).Z();
    }

    m_ResultCurve = GetCurrentFairCurve(theCurve, M, V, D0,D, alpha);
    Iterator(D0, D);
}

void CurveFair::Iterator(Eigen::MatrixXd D0, Eigen::MatrixXd D)
{
    // �ֲ��ƶ�����
    int n = D0.rows();
    Eigen::MatrixXd currentD = D0;
    Standard_Real stepSize = 50;
    Standard_Integer stepCnt = 0;
    bool reached = Standard_False;
    Standard_Real ContorlPointOffset;
    Standard_Real LastControlPointOffset = ControlPointOffsetTol;
    int cnt = 0;
    do
    {
        // ������һ����λ��
        currentD +=  (D - D0) / stepSize;
        cnt++;
        // �̶���β���Ƶ�
        currentD.row(0) << FirstPole.X(), FirstPole.Y(), FirstPole.Z();
        currentD.row(n - 1) << LastPole.X(), LastPole.Y(), LastPole.Z();


        // ����������
        m_ResultCurve = CreateNewBSplineCurve(m_OriginalCurve, currentD);
        ComputeInnerByArcReparam(m_ResultCurve);
        M = ComputeEnergyMatrix(m_OriginalCurve, m_OriginalCurve->Degree(), k);
        Eigen::MatrixXd V(n, n);
        V.setIdentity();
        m_ResultCurve = GetCurrentFairCurve(m_ResultCurve, M, V, D0, D, alpha);
        ContorlPointOffset = GetControlPointsOffset(m_OriginalCurve->Poles(), m_ResultCurve->Poles());
        if (ContorlPointOffset > ControlPointOffsetTol || cnt == stepSize)
        {
            break;
        }

    } while (Standard_True);
}


std::string to_string(Standard_Real value, Standard_Integer precision) 
{
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(precision) << value;
    return oss.str();
}

Handle(Geom_BSplineCurve) CurveFair::GetCurrentFairCurve(const Handle(Geom_BSplineCurve)& theCurve, Eigen::MatrixXd M, Eigen::MatrixXd V, Eigen::MatrixXd D0, Eigen::MatrixXd& D, Standard_Real alpha)
{
    Standard_Real ContorlPointOffset;
    Standard_Real CurveHausdorffDistance;
    Standard_Integer cnt = 0;
    Standard_Integer n = D0.rows();
    Handle(Geom_BSplineCurve) ResultCurve;
    Standard_Real LastControlPointOffset = ControlPointOffsetTol;
    do
    {
        Eigen::MatrixXd temp = (alpha * M + V);
        Eigen::MatrixXd temp_inverse = temp.inverse();

        D = temp_inverse * V * D0;
        D(0, 0) = FirstPole.X();
        D(0, 1) = FirstPole.Y();
        D(0, 2) = FirstPole.Z();
        D(n - 1, 0) = LastPole.X();
        D(n - 1, 1) = LastPole.Y();
        D(n - 1, 2) = LastPole.Z();
        //std::cout << "Alpha" << std::endl << alpha << std::endl;
        //std::cout << "M" << std::endl << format_as(M) << std::endl;
        //std::cout << "V" << std::endl << format_as(V) << std::endl;
        //std::cout << "Alpha * M" << std::endl << format_as(alpha * M) << std::endl;
        //std::cout << "Temp" << std::endl << format_as(temp) << std::endl;
        //std::cout << "Temp_Inverse" << std::endl << format_as(temp_inverse) << std::endl;
        //std::cout << "D0" << std::endl << format_as(Poles) << std::endl;
        //std::cout << "D" << std::endl << format_as(D) << std::endl;


        ResultCurve = CreateNewBSplineCurve(theCurve, D);
        ContorlPointOffset = GetControlPointsOffset(theCurve->Poles(), ResultCurve->Poles());
        CurveHausdorffDistance = GetCurveCurveHausdorffDistance(ResultCurve, theCurve);
        //std::cout << "Cnt:" << cnt << std::endl;
        //std::cout << "Alpha:" << alpha << std::endl;
        //std::cout << "HausdorffDistance:" << CurveHausdorffDistance << std::endl;
        //std::cout << "ControlPointOffset" << ContorlPointOffset << std::endl;
        if (cnt++ > 100)
        {
            Standard_Real Minus = std::abs(ContorlPointOffset - LastControlPointOffset);
            std::cout << "Minus:" << Minus << std::endl;
            if (Minus < 1e-3)
            {
                return ResultCurve;
            }
        }
        if (ContorlPointOffset > ControlPointOffsetTol)
        {
            std::cout << "ContorlPointOffset:" << ContorlPointOffset << std::endl;
            alpha /= 1.5;
        }
        else
        {
            return ResultCurve;
            break;
        }
        LastControlPointOffset = ContorlPointOffset;
    } while (Standard_True);
}
 