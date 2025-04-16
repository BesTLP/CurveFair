#include "curveFair.h"
#include "SurfaceModelingTool.h"
#include "RealCompare.h"
#include "../../OpenCASCADE-7.7.0-vc14-64/opencascade-7.7.0/src/CPnts/CPnts_AbscissaPoint.cxx"
#include <math.hxx>

void UniformCurve(Handle(Geom_BSplineCurve)& curve);
// 定义静态成员变量
std::vector<Handle(Geom_BSplineCurve)> CurveFair::TempResultArray;
std::vector<Handle(Geom_BSplineCurve)> CurveFair::TempInitialArray;
std::vector<Handle(Geom_BSplineCurve)> CurveFair::TempIteratorArray;
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

std::string format_as(gp_Vec Vec)
{
    std::stringstream ss;
    ss << "(" << Vec.X() << ", " << Vec.Y() << ", " << Vec.Z() << ")" << " 模长 : " << Vec.Magnitude();
    return ss.str();
}

std::string to_string(Standard_Real value, Standard_Integer precision)
{
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(precision) << value;
    return oss.str();
}
/*
* 算法: ① 计算出bspline的总长度；② 计算均匀长度step =  bspline / knotNum  - 1；③ 获取均匀长度对应的param；
*/
Standard_Real CurveFair::GetLengthByParam(Standard_Real sParam, Standard_Real eParam, Standard_Real tol)
{
    GeomAdaptor_Curve curveAdapter(outer);
    return GCPnts_AbscissaPoint::Length(curveAdapter, sParam, eParam, tol);
}

Handle(Geom_BSplineCurve) CurveFair::ComputeInnerByArcReparam(const Handle(Geom_BSplineCurve)& theCurve, const Standard_Real tol)
{
    // bspline属性
    Standard_Real firstParam = theCurve->FirstParameter();
    Standard_Real lastParam = theCurve->LastParameter();

    auto bsplineKnots = theCurve->Knots();
    const Standard_Integer LOWER = bsplineKnots.Lower();
    const Standard_Integer UPPER = bsplineKnots.Upper();

    // bspline: 计算弧长
    GeomAdaptor_Curve curveAdapter(theCurve);
    Standard_Real bsplineLen = 0;
    Standard_Real avgLen = 0;
    Standard_Real paramStep = 0;
    try
    {
        bsplineLen = GCPnts_AbscissaPoint::Length(curveAdapter, firstParam, lastParam, tol);
        avgLen = bsplineLen / (paraNum - 1.0);
        paramStep = (lastParam - firstParam) / (paraNum - 1.0);
    }
    catch (const Standard_Failure& e)
    {
        std::cerr << "OpenCASCADE Exception: " << e.GetMessageString() << std::endl;
        return nullptr;
    }
    catch (const std::exception& e)
    {
        std::cerr << "Standard Exception: " << e.what() << std::endl;
        return nullptr;
    }
    catch (...)
    {
        std::cerr << "Unknown Exception occurred!" << std::endl;
        return nullptr;
    }

    // 弧长参数化: 计算新的节点向量
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

    // 为每一个点，添加参数
    fitPoints = reparamPoints;
    /*fitParams = reparamKnots;*/
    inner = GeomAPI_PointsToBSpline
    (
        reparamPoints,
        reparamKnots,
        theCurve->Degree(),
        theCurve->Degree(),
        GeomAbs_C0,
        1e-7
    );
    return inner;
}

Standard_Real CurveFair::f(const Standard_Real theParameter, const Standard_Integer k)
{
    if (k == 0) return inner->Value(theParameter).X();
    gp_Vec DerivateVector = inner->DN(theParameter, k);
    return DerivateVector.X();
}

// To compute the value of a b-spline basic function value （已经验证为正确）
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
 * 计算B样条基函数的k阶导数（已经验证为正确）
 * 基于递归公式: N^{(k)}_{i,p} = p * ((N^{(k-1)}_{i,p-1} / (u_{i+p} - u_i)) - (N^{(k-1)}_{i+1,p-1} / (u_{i+p+1} - u_{i+1})))
 * 可以参考 Nurbs Book 公式(2-9)
 *
 * @param u 参数值
 * @param i 基函数索引
 * @param p 基函数阶数
 * @param k 求导阶数
 * @param Knots 节点向量
 * @return k阶导数值
 */
Standard_Real CurveFair::BasisFunctionDerivative(const Standard_Real u, const Standard_Integer i, const Standard_Integer p,
    const Standard_Integer k, const std::vector<Standard_Real>& Knots)
{
    // 0阶导，直接返回基函数值
    if (k == 0)
    {
        return OneBasicFun(u, i, p, Knots);
    }

    // 当k > p时，导数为0, 当p = 0时，导数为0
    if (k > p || p == 0)
    {
        return 0.0;
    }

    // 计算第一项: (p / (u_{i+p} - u_i)) * N^{(k-1)}_{i,p-1}
    Standard_Real Term1 = 0.0;
    Standard_Real Denominator1 = Knots[i + p] - Knots[i];
    if (std::abs(Denominator1) > std::numeric_limits<Standard_Real>::epsilon())
    {
        Term1 = (p / Denominator1) * BasisFunctionDerivative(u, i, p - 1, k - 1, Knots);
    }

    // 计算第二项: (p / (u_{i+p+1} - u_{i+1})) * N^{(k-1)}_{i+1,p-1}
    Standard_Real Term2 = 0.0;
    Standard_Real Denominator2 = Knots[i + p + 1] - Knots[i + 1];
    if (std::abs(Denominator2) > std::numeric_limits<Standard_Real>::epsilon())
    {
        Term2 = (p / Denominator2) * BasisFunctionDerivative(u, i + 1, p - 1, k - 1, Knots);
    }

    // 根据公式计算k阶导数
    return Term1 - Term2;
}

//Standard_Real CurveFair::GetCurvatureDerivativeSquare(const Standard_Real t, const Standard_Address curF)
//{
//    CurveFair* cf = (CurveFair*)curF;
//    Handle(Geom_BSplineCurve) inner = cf->GetInner();
//    Handle(Geom_BSplineCurve) outer = cf->GetOuter();
//
//    // inner的导数
//    gp_Pnt inPnt;
//    gp_Vec inDerv1, inDerv2, inDerv3;
//    inner->D3(t, inPnt, inDerv1, inDerv2, inDerv3);
//    Standard_Real inParam = inPnt.X();
//    Standard_Real inD1 = inDerv1.X();
//    Standard_Real inD2 = inDerv2.X();
//    Standard_Real inD3 = inDerv3.X();
//
//    // outer的导数
//    gp_Pnt outPnt;
//    gp_Vec outDerv1, outDerv2, outDerv3;
//    outer->D3(inParam, outPnt, outDerv1, outDerv2, outDerv3);
//
//    // D的导数
//    gp_Pnt d0 = outPnt;
//    gp_Vec d1 = outDerv1 * inD1;
//    //std::cout << "单位切向量:" << d1.Magnitude() << std::endl;
//    gp_Vec d2 = outDerv2 * inD1 * inD1 + outDerv1 * inD2;
//    gp_Vec d3 = outDerv3 * pow(inD1, 3) + outDerv2 * (3 * inD1 * inD2) + outDerv1 * inD3;
//    Standard_Real aCurvature = d2.Magnitude() / pow(d1.Magnitude(), 2);
//    Standard_Real  aCurvatureDerivative = d3.Magnitude() / pow(d1.Magnitude(), 3);
//
//    //std::cout << "弧长参数化曲率（修正后）:" << aCurvature << std::endl;
//    //std::cout << "弧长参数化曲率变化率(修正后)" << aCurvatureDerivative << std::endl;
//
//    //// 计算曲率、曲率变化率
//    //Standard_Real aCurvature = d1.Crossed(d2).Magnitude() / pow(d1.Magnitude(), 3);
//    //std::cout << "常规方法曲率:" << aCurvature << std::endl;
//    //gp_Vec dCurvatureDt = (d1.Crossed(d3) * d1.Magnitude() - d1.Crossed(d2) * (3.0 * d1.Dot(d2) / d1.Magnitude())) / pow(d1.Magnitude(), 4);
//    //Standard_Real aCurvatureDerivative = dCurvatureDt.Magnitude() / d1.Magnitude();
//    //std::cout << "常规方法曲率变化率" << aCurvatureDerivative << std::endl;
//    return pow(aCurvatureDerivative, 2);
//}
//
//Standard_Real CurveFair::GetFByGaussIntegral(const Standard_Real tol)
//{
//    CPnts_MyGaussFunction gaussFunction;
//    // POP pout WNT
//    gaussFunction.Init(GetCurvatureDerivativeSquare, this);
//    // FG.Init(f3d,(Standard_Address)&C);
//    // Adaptor3d_Curve adaptor(this->GetOuter());
//    // 使用适配器包装 B-Spline 曲线
//    Handle(GeomAdaptor_Curve) adaptor = new GeomAdaptor_Curve(this->GetOuter());
//
//    std::cout << format_as(this->GetOuter()) << std::endl;
//    std::cout << format_as(this->GetInner()) << std::endl;
//    math_GaussSingleIntegration TheLength(
//        gaussFunction,
//        this->GetOuter()->FirstParameter(),
//        this->GetOuter()->LastParameter(),
//        order(*adaptor),
//        tol
//    );
//
//    return Abs(TheLength.Value());
//}


//struct BasisDerivativeIntegralContext
//{
//    Standard_Integer i;                   // 基函数索引
//    Standard_Integer j;                   // 基函数索引
//    Standard_Integer p;                   // 基函数阶数
//    Standard_Integer k;                   // 导数阶数（对三阶导数为3）
//    std::vector<Standard_Real> Knots;     // 节点向量
//
//    BasisDerivativeIntegralContext
//    (
//        Standard_Integer index_i,
//        Standard_Integer index_j,
//        Standard_Integer degree,
//        Standard_Integer derivOrder,
//        const std::vector<Standard_Real>& knotVector
//    )
//    : i(index_i), j(index_j), p(degree), k(derivOrder), Knots(knotVector) {}
//};
//Standard_Real CurveFair::DerivativeSquaredTwoCallback(const Standard_Real u, const Standard_Address theAddress)
//{
//    // 从地址恢复上下文结构体
//    BasisDerivativeIntegralContext* context =
//        static_cast<BasisDerivativeIntegralContext*>(theAddress);
//
//    Standard_Real fu = f(u);
//    Standard_Real fu_1 = f(u, 1);
//    Standard_Real fu_2 = f(u, 2);
//    Standard_Integer i = context->i;
//    Standard_Integer p = context->p;
//    Standard_Integer DerivateK = context->k;
//    std::vector<Standard_Real> Knots = context->Knots;
//
//
//    Standard_Real Nj3_2 = BasisFunctionDerivative(fu, i, p, DerivateK, Knots);
//    Standard_Real Nj3_1 = BasisFunctionDerivative(fu, i, p, DerivateK, Knots);
//
//    Standard_Real first_1 = Nj3_2 * fu_1 * fu_1;
//    Standard_Real first_2 = Nj3_1 * fu_2;
//    Standard_Real Mi = first_1 + first_2;
//
//    // 返回导数
//    return Mi;
//}
//
//Standard_Real CurveFair::DerivativeSquaredThirdCallback(const Standard_Real u, const Standard_Address theAddress)
//{
//    // 从地址恢复上下文结构体
//    BasisDerivativeIntegralContext* context =
//        static_cast<BasisDerivativeIntegralContext*>(theAddress);
//
//    Standard_Integer i = context->i;
//    Standard_Integer j = context->j;
//
//    Standard_Integer p = context->p;
//    Standard_Integer DerivateK = context->k;
//
//    std::vector<Standard_Real> Knots = context->Knots;
//
//    Standard_Real fu = f(u);
//    Standard_Real fu_1 = f(u, 1);
//    Standard_Real fu_2 = f(u, 2);
//    Standard_Real fu_3 = f(u, 3);
//    Standard_Real Ni3_3 = BasisFunctionDerivative(fu, i, p, DerivateK, Knots);
//    Standard_Real Ni3_2 = BasisFunctionDerivative(fu, i, p, DerivateK, Knots);
//    Standard_Real Ni3_1 = BasisFunctionDerivative(fu, i, p, DerivateK, Knots);
//
//    Standard_Real Nj3_3 = BasisFunctionDerivative(fu, j, p, DerivateK, Knots);
//    Standard_Real Nj3_2 = BasisFunctionDerivative(fu, j, p, DerivateK, Knots);
//    Standard_Real Nj3_1 = BasisFunctionDerivative(fu, j, p, DerivateK, Knots);
//    Standard_Real Mi = Ni3_3 * pow(fu_1, 3) + 3 * Ni3_2 * fu_1 * fu_2 + Ni3_1 * fu_3;
//    Standard_Real Mj = Nj3_3 * pow(fu_1, 3) + 3 * Nj3_2 * fu_1 * fu_2 + Nj3_1 * fu_3;
//    // 返回导数
//    return Mi * Mj;
//}
Standard_Real CurveFair::ComputeDerivativeIntegral
(
    const std::vector<Standard_Real> Knots,
    const Standard_Integer i,
    const Standard_Integer j,
    const Standard_Integer p,
    const Standard_Integer k,
    const Standard_Real tol
)
{
    // 参数有效性验证
    if (i < 0 || i + p + 1 >= Knots.size() || j < 0 || j + p + 1 >= Knots.size() || k > p)
    {
        return 0.0;
    }

    // 确定积分区间 - 基函数的精确非零区间
    // 通常为knots[i] -> kntos[i + p + 1]
    // 计算基函数i和j的非零区间重叠
    Standard_Real lower_i = Knots[i];
    Standard_Real lower_j = Knots[j];
    Standard_Real upper_i = Knots[i + p + 1];
    Standard_Real upper_j = Knots[j + p + 1];

    // 确定有效积分区间
    Standard_Real lowerBound = std::max(lower_i, lower_j);
    Standard_Real upperBound = std::min(upper_i, upper_j);

    //lowerBound = Knots[0];
    //upperBound = Knots[Knots.size() - 1];
    if (upperBound < lowerBound) return 0.0;

    Standard_Real test1 = BasisFunctionDerivative(lowerBound - 0.1, i, p, k, Knots);
    Standard_Real test2 = BasisFunctionDerivative(upperBound + 0.1, i, p, k, Knots);
    Standard_Real test3 = BasisFunctionDerivative(lowerBound - 0.1, j, p, k, Knots);
    Standard_Real test4 = BasisFunctionDerivative(upperBound + 0.1, j, p, k, Knots);
    if ((std::abs(test1) > 1e-10 && std::abs(test3) > 1e-10) ||
          (std::abs(test2) > 1e-10 && std::abs(test4) > 1e-10))
    {
        std::cout << "Error::Integral inteval error!!!" << std::endl;
        std::cout << "test1:" << test1 << std::endl;
        std::cout << "test2:" << test2 << std::endl;
        std::cout << "test3:" << test3 << std::endl;
        std::cout << "test4:" << test4 << std::endl;
    }


    // 检查区间是否有效
    Standard_Real Length = upperBound - lowerBound;
    if (std::abs(upperBound - lowerBound) < 1e-20)
    {
        return 0.0;
    }
    Standard_Integer Step = 100;
    Standard_Real StepParameter = (upperBound - lowerBound) / Step;
    Standard_Real AbsResult = 0;

    //-------------------------------------------------------------------------------
    //                            弧长参数化三阶导数
    //-------------------------------------------------------------------------------

    for (Standard_Real parameter = lowerBound; parameter <= upperBound; parameter += StepParameter)
    {
        if (parameter > upperBound) parameter = upperBound;

        Standard_Real f0 = f(parameter);
        Standard_Real f1 = f(parameter, 1);
        Standard_Real f2 = f(parameter, 2);
        Standard_Real f3 = f(parameter, 3);

        // Ni
        Standard_Real First = BasisFunctionDerivative(f0, i, p, 3, Knots) * pow(f1, 3)
            + 3 * BasisFunctionDerivative(f0, i, p, 2, Knots) * f1 * f2
            + BasisFunctionDerivative(f0, i, p, 1, Knots) * f3;

        // Nj
        Standard_Real Second = BasisFunctionDerivative(f0, j, p, 3, Knots) * pow(f1, 3)
            + 3 * BasisFunctionDerivative(f0, j, p, 2, Knots) * f1 * f2
            + BasisFunctionDerivative(f0, j, p, 1, Knots) * f3;
     
        AbsResult += First * Second;
    }
    return AbsResult * StepParameter;

   ///* -------------------------------------------------------------------------------
   //                             弧长参数化二阶导数
   // -------------------------------------------------------------------------------*/

   // for (Standard_Real parameter = lowerBound; parameter <= upperBound; parameter += StepParameter)
   // {
   //     if (parameter > upperBound) parameter = upperBound;

   //     Standard_Real f0 = f(parameter);
   //     Standard_Real f1 = f(parameter, 1);
   //     Standard_Real f2 = f(parameter, 2);
   //     Standard_Real f3 = f(parameter, 3);

   //     Standard_Real First = BasisFunctionDerivative(f0, i, p, 2, Knots) * pow(f1, 2)
   //         + BasisFunctionDerivative(f0, i, p, 1, Knots) * f2;

   //     Standard_Real Second = BasisFunctionDerivative(f0, j, p, 2, Knots) * pow(f1, 2)
   //         + BasisFunctionDerivative(f0, j, p, 1, Knots) * f2;

   //     AbsResult += First * Second;
   // }
   // return AbsResult * StepParameter;

    
    //-------------------------------------------------------------------------------
    //                             二阶导数
    //-------------------------------------------------------------------------------

    //for (Standard_Real parameter = lowerBound; parameter <= upperBound; parameter += StepParameter)
    //{
    //    if (parameter > upperBound) parameter = upperBound;
    //    Standard_Real First = BasisFunctionDerivative(parameter, i, p, 2, Knots);
    //    Standard_Real Second = BasisFunctionDerivative(parameter, j, p, 2, Knots);

    //    AbsResult += First * Second;
    //}
    //return AbsResult * StepParameter;

    //-------------------------------------------------------------------------------
    //                             三阶导数
    //-------------------------------------------------------------------------------

    //for (Standard_Real parameter = lowerBound; parameter <= upperBound; parameter += StepParameter)
    //{
    //    if (parameter > upperBound) parameter = upperBound;
    //    Standard_Real First = BasisFunctionDerivative(parameter, i, p, 3, Knots);
    //    Standard_Real Second = BasisFunctionDerivative(parameter, j, p, 3, Knots);

    //    AbsResult += First * Second;
    //}
    //return AbsResult * StepParameter;

    //// 创建积分上下文
    //BasisDerivativeIntegralContext context(i, j, p, k, Knots);

    //// 创建高斯积分函数对象并设置回调
    //CPnts_MyGaussFunction gaussFunction;
    //gaussFunction.Init(DerivativeSquaredThirdCallback, &context);
    //// 执行高斯积分
    //math_GaussSingleIntegration integrator(
    //    gaussFunction,
    //    lowerBound,
    //    upperBound,
    //    30,
    //    tol
    //);
    //return Abs(integrator.Value());
}

Standard_Real CurveFair::ComputeDerivativeIntegral
(
    const std::vector<Standard_Real> Knots,
    const Standard_Integer i,
    const Standard_Integer p,
    const Standard_Integer k,
    const Standard_Real tol
)
{
    // 参数有效性验证
    if (i < 0 || i + p + 1 >= Knots.size() || k > p)
    {
        return 0.0;
    }

    // 确定积分区间 - 基函数的精确非零区间
    // 通常为knots[i] -> kntos[i + p + 1]
    // 计算基函数i和j的非零区间重叠
    Standard_Real lowerBound = Knots[i];
    Standard_Real upperBound = Knots[i + p + 1];

    // 检查区间是否有效
    Standard_Real Length = upperBound - lowerBound;
    if (std::abs(upperBound - lowerBound) < 1e-20)
    {
        return 0.0;
    }
    Standard_Integer Step = 100;
    Standard_Real StepParameter = (upperBound - lowerBound) / Step;
    Standard_Real AbsResult = 0;

    //-------------------------------------------------------------------------------
    //                            弧长参数化三阶导数
    //-------------------------------------------------------------------------------

    for (Standard_Real parameter = lowerBound; parameter <= upperBound; parameter += StepParameter)
    {
        if (parameter > upperBound) parameter = upperBound;

        Standard_Real f0 = f(parameter);
        Standard_Real f1 = f(parameter, 1);
        Standard_Real f2 = f(parameter, 2);
        Standard_Real f3 = f(parameter, 3);

        // Ni
        Standard_Real First = BasisFunctionDerivative(f0, i, p, 3, Knots) * pow(f1, 3)
            + 3 * BasisFunctionDerivative(f0, i, p, 2, Knots) * f1 * f2
            + BasisFunctionDerivative(f0, i, p, 1, Knots) * f3;

        AbsResult += First;
    }
    return AbsResult * StepParameter;

    ///* -------------------------------------------------------------------------------
    //                             弧长参数化二阶导数
    // -------------------------------------------------------------------------------*/

    // for (Standard_Real parameter = lowerBound; parameter <= upperBound; parameter += StepParameter)
    // {
    //     if (parameter > upperBound) parameter = upperBound;

    //     Standard_Real f0 = f(parameter);
    //     Standard_Real f1 = f(parameter, 1);
    //     Standard_Real f2 = f(parameter, 2);
    //     Standard_Real f3 = f(parameter, 3);

    //     Standard_Real First = BasisFunctionDerivative(f0, i, p, 2, Knots) * pow(f1, 2)
    //         + BasisFunctionDerivative(f0, i, p, 1, Knots) * f2;

    //     Standard_Real Second = BasisFunctionDerivative(f0, j, p, 2, Knots) * pow(f1, 2)
    //         + BasisFunctionDerivative(f0, j, p, 1, Knots) * f2;

    //     AbsResult += First * Second;
    // }
    // return AbsResult * StepParameter;


     //-------------------------------------------------------------------------------
     //                             二阶导数
     //-------------------------------------------------------------------------------

     //for (Standard_Real parameter = lowerBound; parameter <= upperBound; parameter += StepParameter)
     //{
     //    if (parameter > upperBound) parameter = upperBound;
     //    Standard_Real First = BasisFunctionDerivative(parameter, i, p, 2, Knots);
     //    Standard_Real Second = BasisFunctionDerivative(parameter, j, p, 2, Knots);

     //    AbsResult += First * Second;
     //}
     //return AbsResult * StepParameter;

     //-------------------------------------------------------------------------------
     //                             三阶导数
     //-------------------------------------------------------------------------------

     //for (Standard_Real parameter = lowerBound; parameter <= upperBound; parameter += StepParameter)
     //{
     //    if (parameter > upperBound) parameter = upperBound;
     //    Standard_Real First = BasisFunctionDerivative(parameter, i, p, 3, Knots);
     //    Standard_Real Second = BasisFunctionDerivative(parameter, j, p, 3, Knots);

     //    AbsResult += First * Second;
     //}
     //return AbsResult * StepParameter;

     //// 创建积分上下文
     //BasisDerivativeIntegralContext context(i, j, p, k, Knots);

     //// 创建高斯积分函数对象并设置回调
     //CPnts_MyGaussFunction gaussFunction;
     //gaussFunction.Init(DerivativeSquaredThirdCallback, &context);
     //// 执行高斯积分
     //math_GaussSingleIntegration integrator(
     //    gaussFunction,
     //    lowerBound,
     //    upperBound,
     //    30,
     //    tol
     //);
     //return Abs(integrator.Value());
}

Handle(Geom_BSplineCurve) CurveFair::CreateNewBSplineCurve(const Handle(Geom_BSplineCurve)& originalCurve, const Eigen::MatrixXd& newD)
{
    // 获取原始曲线的 Degree
    Standard_Integer degree = originalCurve->Degree();
    // 获取原始曲线的 节点向量和次数
    TColStd_Array1OfReal knots(1, originalCurve->NbKnots());
    originalCurve->Knots(knots);

    TColStd_Array1OfInteger multiplicities(1, originalCurve->NbKnots());
    originalCurve->Multiplicities(multiplicities);

    // 将 Eigen::Matrix3d 转换为 OCC 控制点数组
    TColgp_Array1OfPnt Poles(1, newD.rows()); // OCC 控制点数组
    for (Standard_Integer i = 0; i < newD.rows(); i++)
    {
        Poles.SetValue(i + 1, gp_Pnt(newD(i, 0), newD(i, 1), newD(i, 2)));
    }
    // 创建新的 B 样条曲线
    Handle(Geom_BSplineCurve) newCurve = new Geom_BSplineCurve(
        Poles, knots, multiplicities, degree);

    UniformCurve(newCurve);

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

        // 检查是否找到投影点
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

Standard_Real CurveFair::f_inverse(const Handle(Geom_BSplineCurve)& theBSplineCurve, Standard_Real t)
{
    Standard_Real Left = theBSplineCurve->FirstParameter();
    Standard_Real Right = theBSplineCurve->LastParameter();
    if (t == Left) return Left;
    if (t == Right) return Right;
    while (true)
    {
        Standard_Real mid = Left + (Right - Left) / 2;
        Standard_Real t_mid = f(mid);
        if (std::abs(t_mid - t) < 1e-10) return mid;
        if (t_mid > t)
        {
            Right = mid;
        }
        else if (t_mid < t)
        {
            Left = mid;
        }
    }
    return 0;
}

Eigen::MatrixXd CurveFair::ComputeEnergyMatrix(
    const Handle(Geom_BSplineCurve)& theBSplineCurve,
    const Standard_Integer p,
    const Standard_Integer k,
    const Standard_Real tol)

{
    const Standard_Integer maxDerivate = k;
    const Standard_Integer aGaussNum = 30;
    Standard_Integer aNum = theBSplineCurve->NbPoles();
    TColStd_Array1OfReal aKnots = theBSplineCurve->Knots();
    Standard_Integer aDeg = theBSplineCurve->Degree();
    TColStd_Array1OfReal aKnotSeq = theBSplineCurve->KnotSequence();
    math_Vector aGaussPnts(1, aGaussNum);
    math_Vector aGausssWeights(1, aGaussNum);
    math::GaussPoints(aGaussNum, aGaussPnts);
    math::GaussWeights(aGaussNum, aGausssWeights);
    std::vector<Standard_Real> KnotVector = ConvertToVector(aKnots);
    std::vector<Standard_Real> KnotSeqVector = ConvertToVector(aKnotSeq);
    Eigen::MatrixXd M(aNum, aNum);
    M.setZero();
    //求解B样条基函数以if*Null的积分
    for (Standard_Integer i = aKnots.Lower(); i < aKnots.Upper(); ++i)
    {
        Standard_Real sStart = f_inverse(theBSplineCurve, aKnots[i]);
        Standard_Real sEnd = f_inverse(theBSplineCurve, aKnots[i + 1]);
        Standard_Real t;
        //将各段的两个高斯积分点相加
        for (Standard_Integer GaussIndex = 0; GaussIndex < aGaussNum; ++GaussIndex)
        {
            t = (sEnd - sStart) * aGaussPnts(aGaussNum - GaussIndex) / 2.0 + (sStart + sEnd) / 2.0;
            Standard_Integer aFirstIndex;
            math_Matrix aBsplineBasis(1, maxDerivate + 1, 1, aDeg + 1);
            BSplCLib::EvalBsplineBasis(maxDerivate, aDeg + 1, aKnotSeq, f(t), aFirstIndex, aBsplineBasis);

            // 计算f(t)及其导数
            Standard_Real f0 = f(t, 0);   // f(t)
            Standard_Real f1 = f(t, 1);   // f'(t)
            Standard_Real f2 = f(t, 2);   // f''(t)
            Standard_Real f3 = f(t, 3);   // f'''(t)
            for (Standard_Integer m = 0; m < aDeg + 1; ++m)
            {
                Standard_Integer globalI = m + aFirstIndex - 1;
                // 原始 k 阶导数 参数 t
                Standard_Real d1Ni = BasisFunctionDerivative(t, globalI, p, 1, KnotSeqVector);
                Standard_Real d2Ni = BasisFunctionDerivative(t, globalI, p, 2, KnotSeqVector);
                Standard_Real d3Ni = BasisFunctionDerivative(t, globalI, p, 3, KnotSeqVector);
                // 弧长参数化 k 阶导数 参数 f(t)
                Standard_Real d1Ni_df1  = BasisFunctionDerivative(f0, globalI, p, 1, KnotSeqVector);
                Standard_Real d2Ni_df2 = BasisFunctionDerivative(f0, globalI, p, 2, KnotSeqVector);
                Standard_Real d3Ni_df3 = BasisFunctionDerivative(f0, globalI, p, 3, KnotSeqVector);

                // 计算弧长参数化链式求导的三阶导数
                Standard_Real D1Ni_Dt1 = d1Ni_df1 * f1;
                Standard_Real D2Ni_Dt2 = d2Ni_df2 * pow(f1, 2) + 3 * d1Ni_df1 * f2;
                Standard_Real D3Ni_Dt3 = d3Ni_df3 * pow(f1, 3) + 3 * d2Ni_df2 * f1 * f2 + d1Ni_df1 * f3;

                for (Standard_Integer n = 0; n < aDeg + 1; ++n)
                {
                    Standard_Integer globalJ = n + aFirstIndex - 1;
                    // 原始 k 阶导数 参数 t
                    Standard_Real d1Nj = BasisFunctionDerivative(t, globalJ, p, 1, KnotSeqVector);
                    Standard_Real d2Nj = BasisFunctionDerivative(t, globalJ, p, 2, KnotSeqVector);
                    Standard_Real d3Nj = BasisFunctionDerivative(t, globalJ, p, 3, KnotSeqVector);
                    // 弧长参数化 k 阶导数 参数 f(t)
                    Standard_Real d1Nj_df1 = BasisFunctionDerivative(f0, globalJ, p, 1, KnotSeqVector);
                    Standard_Real d2Nj_df2 = BasisFunctionDerivative(f0, globalJ, p, 2, KnotSeqVector);
                    Standard_Real d3Nj_df3 = BasisFunctionDerivative(f0, globalJ, p, 3, KnotSeqVector);
                    // 计算弧长参数化链式求导的三阶导数
                    Standard_Real D1Nj_Dt1 = d1Nj_df1 * f1;
                    Standard_Real D2Nj_Dt2 = d2Nj_df2 * pow(f1, 2) + 3 * d1Nj_df1 * f2;
                    Standard_Real D3Nj_Dt3 = d3Nj_df3 * pow(f1, 3) + 3 * d2Nj_df2 * f1 * f2 + d1Nj_df1 * f3;

                    // 弧长参数化三阶导数
                    Standard_Real aElement =
                        D3Ni_Dt3 *
                        D3Nj_Dt3 *
                        aGausssWeights(aGaussNum - GaussIndex) * (sEnd - sStart) / 2.0;
                    M(globalI, globalJ) += aElement;

                }
            }
        }
    }

    //std::cout << M << std::endl;
    //Standard_Real maxVal = M.maxCoeff();
    //if (maxVal > 0)
    //{
    //    M /= maxVal;
    //}
    //std::cout << M << std::endl;
    return M;
}

Standard_Boolean isPositiveDefinite(const Eigen::MatrixXd& A) 
{
    if (A.rows() != A.cols()) return Standard_False;
    Eigen::LLT<Eigen::MatrixXd> llt(A);
    return llt.info() == Eigen::Success;
}

Standard_Boolean isNegativeDefinite(const Eigen::MatrixXd& A) 
{
    // 检查是否对称
    if (!A.isApprox(A.transpose(), 1e-8)) 
    {
        return Standard_False; // 非对称矩阵既不是负定也不是正定
    }

    // 计算特征值
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(A);
    if (solver.info() != Eigen::Success) 
    {
        throw std::runtime_error("Eigen decomposition failed");
    }

    // 检查所有特征值是否 < 0
    return (solver.eigenvalues().array() < 0).all();
}

Handle(Geom_BSplineCurve) CurveFair::GetCurrentFairCurve(
    const Handle(Geom_BSplineCurve)& theCurve,
    Eigen::MatrixXd M,
    Eigen::MatrixXd V,
    Eigen::MatrixXd& D0,
    Eigen::MatrixXd& D,
    Standard_Real alpha)
{
    // 判断是否正定
    std::cout << M.determinant() << std::endl;  // 行列式
    Standard_Boolean M_Positive = isPositiveDefinite(M);
    Standard_Boolean M_Negative = isNegativeDefinite(M);
    Standard_Boolean aM_Positive = isPositiveDefinite(alpha * M + V);
    Standard_Boolean aM_Negative = isNegativeDefinite(alpha * M + V);

    // 控制点的总数
    Standard_Integer n = D0.rows();
    
    // 固定控制点
    Eigen::Vector3d D0_start = D0.row(0);
    Eigen::Vector3d D0_end = D0.row(n - 1);

    // 获取需要优化的控制点原始值
    Eigen::MatrixXd D0_internal = D0.block(1, 0, n - 2, 3);
    Eigen::VectorXd D0_internal_x = D0.col(0).segment(1, n - 2);  // x分量
    Eigen::VectorXd D0_internal_y = D0.col(1).segment(1, n - 2);  // y分量
    Eigen::VectorXd D0_internal_z = D0.col(2).segment(1, n - 2);  // z分量

    // A = α * M_internal + V_internal
    // 获取需要优化的控制点对应的 M 矩阵和 V 矩阵
    Eigen::MatrixXd M_internal = M.block(1, 1, n - 2, n - 2);
    Eigen::MatrixXd V_internal = V.block(1, 1, n - 2, n - 2);

    Standard_Integer ComputeTimes = 0;
    Handle(Geom_BSplineCurve) ResultCurve;
    while (Standard_True)
    {
        ComputeTimes++;
        Eigen::MatrixXd A = alpha * M_internal + V_internal;

        // 分别求解x、y、z分量
        Eigen::VectorXd b_x = Eigen::VectorXd::Zero(n - 2);
        Eigen::VectorXd b_y = Eigen::VectorXd::Zero(n - 2);
        Eigen::VectorXd b_z = Eigen::VectorXd::Zero(n - 2);

        // 减去首尾控制点的影响，分别计算每个分量的b
        for (Standard_Integer k = 0; k < n - 2; k++)
        {
            Standard_Real D0_term_x = M(0, k + 1) * D0_start.x();
            Standard_Real D0_term_y = M(0, k + 1) * D0_start.y();
            Standard_Real D0_term_z = M(0, k + 1) * D0_start.z();
            Standard_Real Dn_term_x = M(n - 1, k + 1) * D0_end.x();
            Standard_Real Dn_term_y = M(n - 1, k + 1) * D0_end.y();
            Standard_Real Dn_term_z = M(n - 1, k + 1) * D0_end.z();

            Standard_Real reg_term_x = V_internal(k, k) * D0_internal_x(k);
            Standard_Real reg_term_y = V_internal(k, k) * D0_internal_y(k);
            Standard_Real reg_term_z = V_internal(k, k) * D0_internal_z(k);

            b_x(k) = -alpha * (D0_term_x + Dn_term_x) + reg_term_x;
            b_y(k) = -alpha * (D0_term_y + Dn_term_y) + reg_term_y;
            b_z(k) = -alpha * (D0_term_z + Dn_term_z) + reg_term_z;
        }

        Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
        Eigen::MatrixXd A_pinv = svd.solve(Eigen::MatrixXd::Identity(A.rows(), A.cols()));
        Eigen::VectorXd D_internal_x = A_pinv * b_x;
        Eigen::VectorXd D_internal_y = A_pinv * b_y;
        Eigen::VectorXd D_internal_z = A_pinv * b_z;

        // 合并结果
        Eigen::MatrixXd D_internal(n - 2, 3);
        D_internal.col(0) = D_internal_x;
        D_internal.col(1) = D_internal_y;
        D_internal.col(2) = D_internal_z;

        // 构建最终控制点矩阵
        D = Eigen::MatrixXd::Zero(n, 3);
        D.row(0) = D0_start;                    // 起始点
        D.block(1, 0, n - 2, 3) = D_internal;   // 内部点
        D.row(n - 1) = D0_end;                  // 终止点

        // 创建新的B样条曲线
        ResultCurve = CreateNewBSplineCurve(m_OriginalCurve, D);

        // 计算新旧曲线的控制点偏差
        Standard_Real ControlPointsOffset = CurveFair::GetControlPointsOffset(m_OriginalCurve->Poles(), ResultCurve->Poles());

        // 计算新曲线的能量
        Standard_Real XEnergy = D.col(0).transpose() * M * D.col(0);
        Standard_Real YEnergy = D.col(1).transpose() * M * D.col(1);
        Standard_Real ZEnergy = D.col(2).transpose() * M * D.col(2);
        std::cout << "XEnergy: " << XEnergy << std::endl;
        std::cout << "YEnergy: " << YEnergy << std::endl;
        std::cout << "ZEnergy: " << ZEnergy << std::endl;
        FairEnergy = alpha *
            (D.col(0).transpose() * M * D.col(0) +
                D.col(1).transpose() * M * D.col(1) +
                D.col(2).transpose() * M * D.col(2)).value();
        std::cout << "-----------" << std::endl;
        std::cout << OriEnergy << std::endl;
        std::cout << FairEnergy << std::endl;
        TempIteratorArray.push_back(ResultCurve);
        if (ControlPointsOffset <= ControlPointOffsetTol || ComputeTimes >= 100)
        {
            break;
        }
        else
        {
            alpha /= 2;
        }

    }
    return ResultCurve;
}

void CurveFair::Iterator(Eigen::MatrixXd D0, Eigen::MatrixXd D)
{
    // 分步移动参数
    Standard_Integer n = D0.rows();
    Standard_Real stepSize = 10;
    Standard_Integer stepCnt = 0;
    Standard_Real ContorlPointOffset;
    Standard_Real LastControlPointOffset = ControlPointOffsetTol;
    while (Standard_True)
    {
        std::cout << D - D0 << std::endl;
        std::cout << "---" << std::endl;
        std::cout << D0 << std::endl;
        std::cout << "---" << std::endl;
        // 计算下一步的位置
        D0 += (D - D0) / 10;
        stepCnt++;

        // 生成新曲线
        Handle(Geom_BSplineCurve) aCurve = CreateNewBSplineCurve(m_OriginalCurve, D0);
        TempInitialArray.push_back(aCurve);

        ComputeInnerByArcReparam(aCurve); // 计算f(t)
        M = ComputeEnergyMatrix(aCurve, m_OriginalCurve->Degree(), k);
        aCurve = GetCurrentFairCurve(aCurve, M, V, D0, D, alpha);
        TempResultArray.push_back(aCurve);

        //ContorlPointOffset = GetControlPointsOffset(m_OriginalCurve->Poles(), aCurve->Poles());
        if (stepCnt == stepSize)
        {
            break;
        }
        else
        {
            m_ResultCurve = aCurve;
        }
    }
}

void CurveFair::Perform(const Handle(Geom_BSplineCurve)& theCurve)
{
    ComputeInnerByArcReparam(theCurve); // 弧长参数化的变换f
    M = ComputeEnergyMatrix(theCurve, theCurve->Degree(), k);

    V.resize(theCurve->NbPoles(), theCurve->NbPoles());
    V.setZero();
    std::vector<Standard_Real> grevilleAbscissae(theCurve->NbPoles());
    for (Standard_Integer i = 0; i < theCurve->NbPoles(); i++)
    {
        Standard_Real sum = 0.0;
        for (Standard_Integer j = 1; j <= theCurve->Degree(); j++)
        {
            Standard_Real aKnotValue = myKnotSeq[i + j];
            aKnotValue = f_inverse(theCurve, aKnotValue);
            sum += aKnotValue;
        }
        grevilleAbscissae[i] = sum / theCurve->Degree();
    }

    std::vector<Standard_Real> curvatureRates(theCurve->NbPoles());
    for (Standard_Integer i = 0; i < theCurve->NbPoles(); i++)
    {
        Standard_Real f0 = f(grevilleAbscissae[i]);
        Standard_Real f1 = f(grevilleAbscissae[i], 1);
        Standard_Real f2 = f(grevilleAbscissae[i], 2);
        Standard_Real f3 = f(grevilleAbscissae[i], 3);

        // Ni
        Standard_Real Curvature = BasisFunctionDerivative(f0, i, theCurve->Degree(), 3, myKnotSeq) * pow(f1, 3)
            + 3 * BasisFunctionDerivative(f0, i, theCurve->Degree(), 2, myKnotSeq) * f1 * f2
            + BasisFunctionDerivative(f0, i, theCurve->Degree(), 1, myKnotSeq) * f3;

        curvatureRates[i] = Curvature * Curvature;
    }

    Standard_Real minRate = *std::min_element(curvatureRates.begin(), curvatureRates.end());
    Standard_Real maxRate = *std::max_element(curvatureRates.begin(), curvatureRates.end());

    for (Standard_Integer i = 0; i < theCurve->NbPoles(); i++)
    {
        // 归一化 curvatureRate 到 [0,1]
        Standard_Real normalized = (curvatureRates[i] - minRate) / (maxRate - minRate + 1e-9);
        // 映射到 [1.0, 0.1]
        Standard_Real weight = 1.0 - 0.5 * normalized; 
        V(i, i) = weight;
    }


    std::cout << V << std::endl;
    TempInitialArray.push_back(theCurve);
    TColgp_Array1OfPnt Poles = theCurve->Poles();
    Standard_Integer n = theCurve->NbPoles();

    Eigen::MatrixXd D0(n, 3); // 原始控制点
    Eigen::MatrixXd D(n, 3); // 新控制点
    for (Standard_Integer i = 0; i < n; i++)
    {
        D0(i, 0) = Poles.Value(i + 1).X();
        D0(i, 1) = Poles.Value(i + 1).Y();
        D0(i, 2) = Poles.Value(i + 1).Z();
    }
    // 初始能量
    OriEnergy = alpha * (D0.col(0).transpose() * M * D0.col(0) + D0.col(1).transpose() * M * D0.col(1) + D0.col(2).transpose() * M * D0.col(2)).value();

    // 获取优化曲线
    m_ResultCurve = GetCurrentFairCurve(theCurve, M, V, D0, D, alpha);
    TempResultArray.push_back(m_ResultCurve);
    Iterator(D0, D);
}