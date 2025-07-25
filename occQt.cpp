/*
*    Copyright (c) 2018 Shing Liu All Rights Reserved.
*
*           File : occQt.cpp
*         Author : Shing Liu(eryar@163.com)
*           Date : 2018-01-08 21:00
*        Version : OpenCASCADE7.2.0 & Qt5.7.1
*
*    Description : Qt main window for OpenCASCADE.
*/

#include "occQt.h"
#include "occView.h"

#include <QToolBar>
#include <QTreeView>
#include <QMessageBox>
#include <QDockWidget>
#include <QProgressBar>

#include <gp_Circ.hxx>
#include <gp_Elips.hxx>
#include <gp_Pln.hxx>

#include <gp_Lin2d.hxx>

#include <Geom_ConicalSurface.hxx>
#include <Geom_ToroidalSurface.hxx>
#include <Geom_CylindricalSurface.hxx>

#include <GCE2d_MakeSegment.hxx>

#include <TopoDS.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TColgp_Array1OfPnt2d.hxx>
#include <Geom_BSplineCurve.hxx>

#include <BRepLib.hxx>

#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>

#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCone.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeTorus.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepPrimAPI_MakeRevol.hxx>

#include <BRepFilletAPI_MakeFillet.hxx>
#include <BRepFilletAPI_MakeChamfer.hxx>

#include <BRepOffsetAPI_MakePipe.hxx>
#include <BRepOffsetAPI_ThruSections.hxx>

#include <iostream>
#include <chrono>
#include <thread>

#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include "D:\OpenCASCADE-7.7.0-vc14-64\opencascade-7.7.0\testPoly\ModelTriangulation.h"
#include <AIS_Shape.hxx>
#include "../../OpenCASCADE-7.7.0-vc14-64/qt5.11.2-vc14-64/include/QtWidgets/qfiledialog.h"
#include "../../OpenCASCADE-7.7.0-vc14-64/opencascade-7.7.0/src/STEPControl/STEPControl_Reader.hxx"
#include "../../OpenCASCADE-7.7.0-vc14-64/opencascade-7.7.0/src/BRepAlgoAPI/BRepAlgoAPI_Section.hxx"
#include "D:\OpenCASCADE-7.7.0-vc14-64\opencascade-7.7.0\testPoly\ModelImporter.h"
#include "../../OpenCASCADE-7.7.0-vc14-64/qt5.11.2-vc14-64/include/QtWidgets/qprogressdialog.h"
#include "D:\OpenCASCADE-7.7.0-vc14-64\opencascade-7.7.0\testPoly\RandomExport.h"
#include "../../OpenCASCADE-7.7.0-vc14-64/qt5.11.2-vc14-64/include/QtWidgets/qinputdialog.h"
#include <BRepTools.hxx>
#include <Geom_Hyperbola.hxx>
#include <Geom_TrimmedCurve.hxx>
#include "../../OpenCASCADE-7.7.0-vc14-64/opencascade-7.7.0/src/BRepClass3d/BRepClass3d_SolidClassifier.hxx"
#include "../../OpenCASCADE-7.7.0-vc14-64/opencascade-7.7.0/src/BRepExtrema/BRepExtrema_DistShapeShape.hxx"
#include "../../OpenCASCADE-7.7.0-vc14-64/qt5.11.2-vc14-64/include/QtWidgets/qboxlayout.h"
#include "D:\OpenCASCADE-7.7.0-vc14-64\opencascade-7.7.0\testPoly\MakeShape.h"
#include "../../OpenCASCADE-7.7.0-vc14-64/opencascade-7.7.0/testPoly/ModelExporter.h"
#include <future>
#include <mutex>
#include "../../OpenCASCADE-7.7.0-vc14-64/opencascade-7.7.0/src/GeomLib/GeomLib.hxx"
#include <GeomAPI_PointsToBSpline.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Compound.hxx>
#include <BRep_Builder.hxx>
#include <BRepTools.hxx>
#include <BRepTools_ReShape.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <IGESControl_Reader.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Vertex.hxx>
#include <gp_Pnt.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <vector>
#include <limits>
#include <iostream>
#include "../src/IGESCAFControl/IGESCAFControl_Reader.hxx"
#include "../src/BRepBuilderAPI/BRepBuilderAPI_MakeVertex.hxx"
#include "../inc/BRepAdaptor_Curve.hxx"
#include "../src/GCPnts/GCPnts_TangentialDeflection.hxx"
#include "../inc/BRepAdaptor_Surface.hxx"
#include "BRepAdaptor_Surface.hxx"
#include "../inc/STEPControl_Reader.hxx"
#include "Geom_TrimmedCurve.hxx"
#include "../src/GeomLib/GeomLib.hxx"
#include <Standard.hxx>
#include <Standard_DefineAlloc.hxx>

#include <math_DoubleTab.hxx>
#include <math_Vector.hxx>
#include <Standard_OStream.hxx>
#include <math_Matrix.hxx>
#include "../src/Convert/Convert_ParameterisationType.hxx"
#include "../src/TColgp/TColgp_Array1OfXYZ.hxx"
#include "../src/GeomConvert/GeomConvert_CompCurveToBSplineCurve.hxx"
#include "../src/GeomLib/GeomLib.cxx"
#include <TColgp_Array1OfXYZ.hxx>
#include "../../OpenCASCADE-7.7.0-vc14-64/opencascade-7.7.0/src/GeomAPI/GeomAPI_ExtremaCurveCurve.hxx"
#include "../../OpenCASCADE-7.7.0-vc14-64/opencascade-7.7.0/src/GeomAPI/GeomAPI_ProjectPointOnCurve.hxx"
#include "SurfaceModelingTool.h"
#include "../../OpenCASCADE-7.7.0-vc14-64/opencascade-7.7.0/src/GCPnts/GCPnts_QuasiUniformAbscissa.hxx"
#include <GCPnts_AbscissaPoint.hxx>

// Thread
#include <thread>
#include <chrono>
#include <atomic>
#include "../../OpenCASCADE-7.7.0-vc14-64/opencascade-7.7.0/src/GeomFill/GeomFill_Coons.hxx"
#include "BRepBuilderAPI_MakeFace.hxx"
#include <STEPControl_Writer.hxx>
#include <GeomFill_BSplineCurves.hxx>
#include <AIS_TextLabel.hxx>
#include <GordenSurface.h>
#include "LogPrint.h"
#include <CurveFair.h>
#include "D:\OpenCASCADE-7.7.0-vc14-64\opencascade-7.7.0\testPoly\MakeShape.h"
#include <random>

template <typename T>
void occQt::Visualize(const T& object, const Quantity_Color& color)
{
    if (isVisualize == Standard_False) return;
    Handle(AIS_InteractiveContext) context = myOccView->getContext();  // 直接访问 myOccView
    try
    {
        Handle(AIS_Shape) aisShape;

        // 根据对象类型创建对应的 AIS_Shape
        if constexpr (std::is_same<T, gp_Pnt>::value)
        {
            TopoDS_Vertex ver = BRepBuilderAPI_MakeVertex(object);
            aisShape = new AIS_Shape(ver);
        }
        else if constexpr (std::is_same<T, Handle(Geom_BSplineCurve)>::value)
        {
            TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(object);
            aisShape = new AIS_Shape(edge);
        }
        else if constexpr (std::is_same<T, Handle(Geom_BSplineSurface)>::value)
        {
            Handle(Geom_Surface) genericSurface = Handle(Geom_Surface)::DownCast(object);
            if (genericSurface.IsNull())
            {
                return;  // 跳过无效的面
            }
            TopoDS_Face face = BRepBuilderAPI_MakeFace(genericSurface, Precision::Confusion());
            aisShape = new AIS_Shape(face);
        }
        else if constexpr (std::is_same<T, TopoDS_Shape>::value)
        {
            aisShape = new AIS_Shape(object);
        }
        else if constexpr (std::is_same<T, TopoDS_Edge>::value)
        {
            aisShape = new AIS_Shape(object);
        }
        else if constexpr (std::is_same<T, std::pair<gp_Pnt, Standard_Real>>::value)
        {
            const gp_Pnt& point = object.first;
            Standard_Real value = object.second;

            // 格式化Standard_Real，保留两位小数
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(2) << value;
            TCollection_AsciiString textString(oss.str().c_str());

            // 创建文本对象
            Handle(AIS_TextLabel) aisTextLabel = new AIS_TextLabel();
            aisTextLabel->SetText(textString);
            aisTextLabel->SetPosition(point);
            aisTextLabel->SetColor(color);
            // 显示文本
            context->Display(aisTextLabel, Standard_True);
            TopoDS_Vertex ver = BRepBuilderAPI_MakeVertex(point);

            // 显示点
            aisShape = new AIS_Shape(ver);
            aisShape->SetColor(color); // 设置颜色
            context->Display(aisShape, Standard_True);
            return;  // 跳过后续的 AIS_Shape 处理
        }
        else if constexpr (std::is_same<T, gp_Pln>::value)
        {
            // 定义平面的范围
            Standard_Real uMin = -500.0, uMax = 500.0, vMin = -500.0, vMax = 500.0;

            // 将 gp_Pln 转换为 Geom_Plane
            Handle(Geom_Plane) geomPlane = new Geom_Plane(object);

            // 创建 TopoDS_Face
            TopoDS_Face face = BRepBuilderAPI_MakeFace(geomPlane, uMin, uMax, vMin, vMax, Precision::Confusion());
            aisShape = new AIS_Shape(face);
        }

        aisShape->SetColor(color); // 设置颜色
        context->Display(aisShape, Standard_True);
    }
    catch (Standard_Failure& e)
    {
        QMessageBox::critical(nullptr, "Error", QString("Failed to display object: %1").arg(e.GetMessageString()));
        return;
    }

    // 调整视角，使所有对象适应视图
    if (myOccView)
    {
        myOccView->fitAll();
    }
}

template <typename T>
void occQt::Visualize(const std::vector<T>& objects, const Quantity_Color& color)
{
    if (isVisualize == Standard_False) return;
    Handle(AIS_InteractiveContext) context = myOccView->getContext();  // 直接访问 myOccView
    for (const auto& obj : objects)
    {
        try
        {
            Handle(AIS_Shape) aisShape;

            // 根据对象类型创建对应的 AIS_Shape
            if constexpr (std::is_same<T, gp_Pnt>::value)
            {
                TopoDS_Vertex ver = BRepBuilderAPI_MakeVertex(obj);
                aisShape = new AIS_Shape(ver);
            }
            else if constexpr (std::is_same<T, Handle(Geom_BSplineCurve)>::value)
            {
                TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(obj);
                aisShape = new AIS_Shape(edge);
            }
            else if constexpr (std::is_same<T, Handle(Geom_BSplineSurface)>::value)
            {
                Handle(Geom_Surface) genericSurface = Handle(Geom_Surface)::DownCast(obj);
                if (genericSurface.IsNull())
                {
                    continue;
                }
                TopoDS_Face face = BRepBuilderAPI_MakeFace(genericSurface, Precision::Confusion());
                aisShape = new AIS_Shape(face);
            }
            else if constexpr (std::is_same<T, TopoDS_Shape>::value)
            {
                aisShape = new AIS_Shape(obj);
            }
            else if constexpr (std::is_same<T, TopoDS_Edge>::value)
            {
                aisShape = new AIS_Shape(obj);
            }
            else if constexpr (std::is_same<T, std::pair<gp_Pnt, Standard_Real>>::value)
            {
                const gp_Pnt& point = obj.first;
                Standard_Real value = obj.second;

                // 格式化Standard_Real，保留两位小数
                std::ostringstream oss;
                oss << std::fixed << std::setprecision(2) << value;
                TCollection_AsciiString textString(oss.str().c_str());

                Handle(AIS_TextLabel) aisTextLabel = new AIS_TextLabel();
                aisTextLabel->SetText(textString);
                aisTextLabel->SetPosition(point);
                aisTextLabel->SetColor(color);
                // 显示文本
                context->Display(aisTextLabel, Standard_True);
                TopoDS_Vertex ver = BRepBuilderAPI_MakeVertex(point);

                // 显示点
                aisShape = new AIS_Shape(ver);
                aisShape->SetColor(color); // 设置颜色
                context->Display(aisShape, Standard_True);
                continue;  // 跳过后续的 AIS_Shape 处理
            }
            else if constexpr (std::is_same<T, gp_Pln>::value)
            {
                // 定义平面的范围
                Standard_Real uMin = -500.0, uMax = 500.0, vMin = -500.0, vMax = 500.0;

                // 将 gp_Pln 转换为 Geom_Plane
                Handle(Geom_Plane) geomPlane = new Geom_Plane(obj);

                // 创建 TopoDS_Face
                TopoDS_Face face = BRepBuilderAPI_MakeFace(geomPlane, uMin, uMax, vMin, vMax, Precision::Confusion());
                aisShape = new AIS_Shape(face);
            }

            aisShape->SetColor(color); // 设置颜色
            context->Display(aisShape, Standard_True);
        }
        catch (Standard_Failure& e)
        {
            QMessageBox::critical(nullptr, "Error", QString("Failed to display object: %1").arg(e.GetMessageString()));
        }
    }

    // 调整视角，使所有对象适应视图
    if (myOccView)
    {
        myOccView->fitAll();
    }
}


void ExportBSplineSurface(const Handle(Geom_BSplineSurface)& bsplineSurface, const std::string& filename)
{
    // 获取曲面的参数范围
    Standard_Real uMin, uMax, vMin, vMax;
    bsplineSurface->Bounds(uMin, uMax, vMin, vMax);

    // 使用曲面和参数范围创建面
    TopoDS_Face face = BRepBuilderAPI_MakeFace(bsplineSurface, uMin, uMax, vMin, vMax, 1e-7);

    // 检查面是否有效
    if (face.IsNull())
    {
        std::cerr << "面创建失败！" << std::endl;
        return;
    }

    // 将文件名转换为小写以进行不区分大小写的比较
    std::string filename_lower = filename;
    std::transform(filename_lower.begin(), filename_lower.end(), filename_lower.begin(), ::tolower);

    // 根据文件扩展名选择输出格式
    if (filename_lower.size() >= 5 && filename_lower.substr(filename_lower.size() - 5) == ".brep")
    {
        // 将面保存到 BREP 文件
        if (BRepTools::Write(face, filename.c_str()))
        {
            std::cout << "成功导出到 BREP 文件: " << filename << std::endl;
        }
        else
        {
            std::cerr << "导出 BREP 文件失败！" << std::endl;
        }
    }
    else if (filename_lower.size() >= 5 && filename_lower.substr(filename_lower.size() - 5) == ".step")
    {
        // 将面保存到 STEP 文件
        STEPControl_Writer writer;
        IFSelect_ReturnStatus status = writer.Transfer(face, STEPControl_AsIs);
        if (status == IFSelect_RetDone)
        {
            status = writer.Write(filename.c_str());
            if (status == IFSelect_RetDone)
            {
                std::cout << "成功导出到 STEP 文件: " << filename << std::endl;
            }
            else
            {
                std::cerr << "导出 STEP 文件失败！" << std::endl;
            }
        }
        else
        {
            std::cerr << "面转换为 STEP 格式失败！" << std::endl;
        }
    }
    else
    {
        std::cerr << "不支持的文件扩展名，请使用 .brep 或 .step。" << std::endl;
    }
}


void PrintMatrix(const math_Matrix& Mat) {
    for (Standard_Integer i = Mat.LowerRow(); i <= Mat.UpperRow(); ++i) {
        for (Standard_Integer j = Mat.LowerCol(); j <= Mat.UpperCol(); ++j) {
            std::cout << Mat(i, j) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "----------------------------" << std::endl;
}

void PrintArray1OfXYZ(const TColgp_Array1OfXYZ& array)
{
    for (Standard_Integer i = array.Lower(); i <= array.Upper(); ++i)
    {
        const gp_XYZ& point = array.Value(i);
        std::cout << "Point " << i << ": ("
            << point.X() << ", "
            << point.Y() << ", "
            << point.Z() << ")" << std::endl;
    }
    std::cout << "---------------------" << std::endl;
}
void PrintArray1OfPnt(const TColgp_Array1OfPnt& array)
{
    for (Standard_Integer i = array.Lower(); i <= array.Upper(); ++i)
    {
        const gp_Pnt& point = array.Value(i);
        std::cout << "Point " << i << ": ("
            << point.X() << ", "
            << point.Y() << ", "
            << point.Z() << ")" << std::endl;
    }
    std::cout << "---------------------" << std::endl;
}
void M_ExtendCurveToPoint(Handle(Geom_BoundedCurve)& Curve,
    const gp_Pnt& Point,
    const Standard_Integer Continuity,
    const Standard_Boolean After)
{
    if (Continuity < 1 || Continuity > 3) return;
    Standard_Integer size = Continuity + 2;
    Standard_Real Ubord, Tol = 1.e-6;
    math_Matrix  MatCoefs(1, size, 1, size);
    Standard_Real Lambda, L1;
    Standard_Integer ii, jj;
    gp_Vec d1, d2, d3;
    gp_Pnt p0;
    // il faut Convertir l'entree (en preservant si possible le parametrage)
    GeomConvert_CompCurveToBSplineCurve Concat(Curve, Convert_QuasiAngular);

    // Les contraintes de constructions
    TColgp_Array1OfXYZ Cont(1, size);
    if (After) {
        Ubord = Curve->LastParameter();

    }
    else {
        Ubord = Curve->FirstParameter();
    }
    PLib::HermiteCoefficients(0, 1,           // Les Bornes
        Continuity, 0,  // Les Ordres de contraintes
        MatCoefs);
    PrintMatrix(MatCoefs);

    Curve->D3(Ubord, p0, d1, d2, d3);
    if (!After) { // Inversion du parametrage
        d1 *= -1;
        d3 *= -1;
    }

    L1 = p0.Distance(Point);
    if (L1 > Tol) {
        // Lambda est le ratio qu'il faut appliquer a la derive de la courbe
        // pour obtenir la derive du prolongement (fixe arbitrairement a la
        // longueur du segment bout de la courbe - point cible.
        // On essai d'avoir sur le prolongement la vitesse moyenne que l'on
        // a sur la courbe.
        gp_Vec daux;
        gp_Pnt pp;
        Standard_Real f = Curve->FirstParameter(), t, dt, norm;
        dt = (Curve->LastParameter() - f) / 9;
        norm = d1.Magnitude();
        for (ii = 1, t = f + dt; ii <= 8; ii++, t += dt) {
            Curve->D1(t, pp, daux);
            norm += daux.Magnitude();
        }
        norm /= 9;
        dt = d1.Magnitude() / norm;
        if ((dt < 1.5) && (dt > 0.75)) { // Le bord est dans la moyenne on le garde
            Lambda = ((Standard_Real)1) / Max(d1.Magnitude() / L1, Tol);
        }
        else {
            Lambda = ((Standard_Real)1) / Max(norm / L1, Tol);
        }
    }
    else {
        return; // Pas d'extension
    }

    // Optimisation du Lambda
    math_Matrix Cons(1, 3, 1, size);
    Cons(1, 1) = p0.X();  Cons(2, 1) = p0.Y(); Cons(3, 1) = p0.Z();
    Cons(1, 2) = d1.X();  Cons(2, 2) = d1.Y(); Cons(3, 2) = d1.Z();
    Cons(1, size) = Point.X();  Cons(2, size) = Point.Y(); Cons(3, size) = Point.Z();
    if (Continuity >= 2) {
        Cons(1, 3) = d2.X();  Cons(2, 3) = d2.Y(); Cons(3, 3) = d2.Z();
    }
    if (Continuity >= 3) {
        Cons(1, 4) = d3.X();  Cons(2, 4) = d3.Y(); Cons(3, 4) = d3.Z();
    }
    PrintMatrix(MatCoefs);
    // Construction dans la Base Polynomiale
    Cont(1) = p0.XYZ();
    Cont(2) = d1.XYZ() * Lambda;
    if (Continuity >= 2) Cont(3) = d2.XYZ() * Pow(Lambda, 2);
    if (Continuity >= 3) Cont(4) = d3.XYZ() * Pow(Lambda, 3);
    Cont(size) = Point.XYZ();
    PrintArray1OfXYZ(Cont);

    TColgp_Array1OfPnt ExtrapPoles(1, size);
    TColgp_Array1OfPnt ExtraCoeffs(1, size);

    gp_Pnt PNull(0., 0., 0.);
    ExtraCoeffs.Init(PNull);
    for (ii = 1; ii <= size; ii++) {
        for (jj = 1; jj <= size; jj++) {
            ExtraCoeffs(jj).ChangeCoord() += MatCoefs(ii, jj) * Cont(ii);
        }
    }
    PrintArray1OfPnt(ExtraCoeffs);
    // Convertion Dans la Base de Bernstein
    PLib::CoefficientsPoles(ExtraCoeffs, PLib::NoWeights(),
        ExtrapPoles, PLib::NoWeights());


    PrintArray1OfPnt(ExtrapPoles);
    Handle(Geom_BezierCurve) Bezier = new (Geom_BezierCurve) (ExtrapPoles);

    Standard_Real dist = ExtrapPoles(1).Distance(p0);
    Standard_Boolean Ok;
    Tol += dist;

    // Concatenation
    Ok = Concat.Add(Bezier, Tol, After);
    if (!Ok) throw Standard_ConstructionError("ExtendCurveToPoint");

    Curve = Concat.BSplineCurve();
}
gp_Ax2 getCoordinateSystemFromUserInput(bool ok)
{
    Standard_Real xOrigin = QInputDialog::getDouble(nullptr, "Input", "Enter X coordinate for origin:", 0, -1000, 1000, 8, &ok);
    Standard_Real yOrigin = QInputDialog::getDouble(nullptr, "Input", "Enter Y coordinate for origin:", 0, -1000, 1000, 8, &ok);
    Standard_Real zOrigin = QInputDialog::getDouble(nullptr, "Input", "Enter Z coordinate for origin:", 0, -1000, 1000, 8, &ok);

    Standard_Real zDirX = QInputDialog::getDouble(nullptr, "Input", "Enter X coordinate for Z direction:", 0, -1000, 1000, 8, &ok);
    Standard_Real zDirY = QInputDialog::getDouble(nullptr, "Input", "Enter Y coordinate for Z direction:", 0, -1000, 1000, 8, &ok);
    Standard_Real zDirZ = QInputDialog::getDouble(nullptr, "Input", "Enter Z coordinate for Z direction:", 0, -1000, 1000, 8, &ok);

    Standard_Real xDirX = QInputDialog::getDouble(nullptr, "Input", "Enter X coordinate for X direction:", 0, -1000, 1000, 8, &ok);
    Standard_Real xDirY = QInputDialog::getDouble(nullptr, "Input", "Enter Y coordinate for X direction:", 0, -1000, 1000, 8, &ok);
    Standard_Real xDirZ = QInputDialog::getDouble(nullptr, "Input", "Enter Z coordinate for X direction:", 0, -1000, 1000, 8, &ok);

    // Create gp_Ax2 with user-provided values
    gp_Ax2 coordinateSystemAx2 = gp_Ax2(gp_Pnt(xOrigin, yOrigin, zOrigin), gp_Dir(zDirX, zDirY, zDirZ), gp_Dir(xDirX, xDirY, xDirZ));

    return coordinateSystemAx2;
}
// Function to create different types of curves
TopoDS_Shape createCurve(const QString& curveType)
{
    bool ok;
    if (curveType.toLower() == "circle")
    {
        // Create a circle here
        // ...
        // Create gp_Ax2 with user-provided values
        gp_Ax2 circleAx2 = getCoordinateSystemFromUserInput(ok);
        // Get user input for major and minor radii
        Standard_Real radius = QInputDialog::getDouble(nullptr, "Input", "Enter Radius:", 0, -1000, 1000, 8, &ok);

        gp_Circ cir(circleAx2, radius);
        BRepBuilderAPI_MakeEdge edgeMaker(cir);
        if (edgeMaker.IsDone())
        {
            return edgeMaker.Shape();
        }
    }
    else if (curveType.toLower() == "ellipse")
    {
        // Create gp_Ax2 with user-provided values
        gp_Ax2 ellipseAx2 = getCoordinateSystemFromUserInput(ok);
        // Get user input for major and minor radii
        Standard_Real majorRadius = QInputDialog::getDouble(nullptr, "Input", "Enter Major Radius:", 0, -1000, 1000, 8, &ok);
        Standard_Real minorRadius = QInputDialog::getDouble(nullptr, "Input", "Enter Minor Radius:", 0, -1000, 1000, 8, &ok);



        gp_Elips ellipse(ellipseAx2, majorRadius, minorRadius);
        BRepBuilderAPI_MakeEdge edgeMaker(ellipse);
        if (edgeMaker.IsDone())
        {
            return edgeMaker.Shape();
        }
    }
    else if (curveType.toLower() == "parabola")
    {
        // Create a parabola here
        // ...
                // Get user input for origin coordinates
        // Create gp_Ax2 with user-provided values
        gp_Ax2 parabolaAx2 = getCoordinateSystemFromUserInput(ok);
        // Get user input for major and minor radii
        Standard_Real theFocal = QInputDialog::getDouble(nullptr, "Input", "Enter the Focal:", 0, -1000, 1000, 8, &ok);

        gp_Parab parabola(parabolaAx2, theFocal);
        BRepBuilderAPI_MakeEdge edgeMaker(parabola);
        if (edgeMaker.IsDone())
        {
            return edgeMaker.Shape();
        }

    }
    else if (curveType.toLower() == "hyperbola")
    {
        // Get user input for origin coordinates
        gp_Ax2 hyperbolaAx2 = getCoordinateSystemFromUserInput(ok);
        // Get user input for major and minor radii
        Standard_Real majorRadius = QInputDialog::getDouble(nullptr, "Input", "Enter Major Radius:", 0, -1000, 1000, 8, &ok);
        Standard_Real minorRadius = QInputDialog::getDouble(nullptr, "Input", "Enter Minor Radius:", 0, -1000, 1000, 8, &ok);

        // 在正半轴上创建双曲线
        gp_Hypr hyperbolaPositive(hyperbolaAx2, majorRadius, minorRadius);
        BRepBuilderAPI_MakeEdge edgeMakerPositive(hyperbolaPositive);

        // 沿着 y 轴进行 90° 的旋转，得到负半轴上的双曲线
        gp_Trsf rotationTrsf;
        gp_Ax1 axis(hyperbolaAx2.Location(), hyperbolaAx2.YDirection());

        rotationTrsf.SetRotation(axis, M_PI);  // M_PI 是圆周率
        gp_Hypr hyperbolaNegative = hyperbolaPositive.Transformed(rotationTrsf);
        BRepBuilderAPI_MakeEdge edgeMakerNegative(hyperbolaNegative);

        // 将正负半轴上的双曲线放入一个复合形状中
        TopoDS_Compound compound;
        BRep_Builder compoundBuilder;
        compoundBuilder.MakeCompound(compound);
        compoundBuilder.Add(compound, edgeMakerPositive.Shape());
        compoundBuilder.Add(compound, edgeMakerNegative.Shape());

        return compound;
    }
    else if (curveType.toLower() == "line")
    {
        // Create a line here
        // ...
        // Get user input for origin coordinates
        Standard_Real xStart = QInputDialog::getDouble(nullptr, "Input", "Enter X coordinate for Start Point:", 0, -1000, 1000, 8, &ok);
        Standard_Real yStart = QInputDialog::getDouble(nullptr, "Input", "Enter Y coordinate for Start Point:", 0, -1000, 1000, 8, &ok);
        Standard_Real zStart = QInputDialog::getDouble(nullptr, "Input", "Enter Z coordinate for Start Point:", 0, -1000, 1000, 8, &ok);

        // Get user input for Z direction coordinates
        Standard_Real DirX = QInputDialog::getDouble(nullptr, "Input", "Enter X coordinate for direction:", 0, -1000, 1000, 8, &ok);
        Standard_Real DirY = QInputDialog::getDouble(nullptr, "Input", "Enter Y coordinate for direction:", 0, -1000, 1000, 8, &ok);
        Standard_Real DirZ = QInputDialog::getDouble(nullptr, "Input", "Enter Z coordinate for direction:", 0, -1000, 1000, 8, &ok);

        gp_Pnt startPoint(xStart, yStart, zStart);
        gp_Dir dir(DirX, DirY, DirZ);
        gp_Lin line(startPoint, dir);
        BRepBuilderAPI_MakeEdge edgeMaker(line);
        if (edgeMaker.IsDone())
        {
            return edgeMaker.Shape();
        }
    }

    // Return an empty shape if the curve type is not recognized
    return TopoDS_Shape();
}
TopoDS_Shape createSurface(const QString& surfaceType)
{
    bool ok;
    if (surfaceType.toLower() == "plane")
    {
        gp_Ax2 planeAxis = getCoordinateSystemFromUserInput(ok);
        Standard_Real length = QInputDialog::getDouble(nullptr, "Input", "Enter length:", 0, -1000, 1000, 8, &ok);
        Standard_Real width = QInputDialog::getDouble(nullptr, "Input", "Enter width:", 0, -1000, 1000, 8, &ok);

        return MakeShape::MakePlane(planeAxis.Location(), planeAxis.XDirection() ^ planeAxis.YDirection(), planeAxis.XDirection(), length, width);

    }
    else if (surfaceType.toLower() == "cylinder")
    {
        // Create gp_Ax2 with user-provided values
        gp_Ax2 cylinderAxis = getCoordinateSystemFromUserInput(ok);
        // Get user input for major and minor radii
        Standard_Real radius = QInputDialog::getDouble(nullptr, "Input", "Enter Radius:", 0, -1000, 1000, 8, &ok);
        Standard_Real height = QInputDialog::getDouble(nullptr, "Input", "Enter Height:", 0, -1000, 1000, 8, &ok);

        return MakeShape::MakeCylinder(cylinderAxis.Location(), cylinderAxis.XDirection() ^ cylinderAxis.YDirection(), cylinderAxis.XDirection(), radius, height);
    }
    else if (surfaceType.toLower() == "sphere")
    {
        gp_Ax2 sphereAxis = getCoordinateSystemFromUserInput(ok);
        // Get user input for major and minor radii
        Standard_Real radius = QInputDialog::getDouble(nullptr, "Input", "Enter Radius:", 0, -1000, 1000, 8, &ok);

        return MakeShape::MakeSphere(sphereAxis.Location(), sphereAxis.XDirection() ^ sphereAxis.YDirection(), sphereAxis.XDirection(), radius);

    }
    else if (surfaceType.toLower() == "cone")
    {
        // Create gp_Ax2 with user-provided values
        gp_Ax2 coneAxis = getCoordinateSystemFromUserInput(ok);
        // Get user input for major and minor radii
        Standard_Real radius = QInputDialog::getDouble(nullptr, "Input", "Enter Radius:", 0, -1000, 1000, 8, &ok);
        Standard_Real height = QInputDialog::getDouble(nullptr, "Input", "Enter Height:", 0, -1000, 1000, 8, &ok);

        return MakeShape::MakeCone(coneAxis.Location(), coneAxis.XDirection() ^ coneAxis.YDirection(), coneAxis.XDirection(), radius, height);
    }
    else if (surfaceType.toLower() == "torus")
    {
        // Create gp_Ax2 with user-provided values
        gp_Ax2 torusAxis = getCoordinateSystemFromUserInput(ok);
        // Get user input for major and minor radii
        Standard_Real majorRadius = QInputDialog::getDouble(nullptr, "Input", "Enter Major Radius:", 0, -1000, 1000, 8, &ok);
        Standard_Real minorRadius = QInputDialog::getDouble(nullptr, "Input", "Enter Minor Radius:", 0, -1000, 1000, 8, &ok);

        return MakeShape::MakeTorus(torusAxis.Location(), torusAxis.XDirection() ^ torusAxis.YDirection(), torusAxis.XDirection(), majorRadius, minorRadius);
    }

    // Return an empty shape if the curve type is not recognized
    return TopoDS_Shape();
}
occQt::occQt(QWidget* parent)
    : QMainWindow(parent)
{
    ui.setupUi(this);

    createActions();
    createMenus();
    createToolBars();

    myOccView = new OccView(this);
    setCentralWidget(myOccView);

}

occQt::~occQt()
{

}

void occQt::createActions(void)
{
    // File
    connect(ui.actionExit, SIGNAL(triggered()), this, SLOT(close()));


    // View
    connect(ui.actionZoom, SIGNAL(triggered()), myOccView, SLOT(zoom()));
    connect(ui.actionPan, SIGNAL(triggered()), myOccView, SLOT(pan()));
    connect(ui.actionRotate, SIGNAL(triggered()), myOccView, SLOT(rotate()));

    connect(ui.actionReset, SIGNAL(triggered()), myOccView, SLOT(reset()));
    connect(ui.actionFitAll, SIGNAL(triggered()), myOccView, SLOT(fitAll()));

    // Primitive
    connect(ui.actionBox, SIGNAL(triggered()), this, SLOT(importFile()));
    connect(ui.actionCone, SIGNAL(triggered()), this, SLOT(Triangulation()));
    connect(ui.actionSphere, SIGNAL(triggered()), this, SLOT(TriangulationIntersection()));
    connect(ui.actionCylinder, SIGNAL(triggered()), this, SLOT(ClearDisplay()));
    connect(ui.actionTorus, SIGNAL(triggered()), this, SLOT(PrintInfo()));

    // Modeling
    connect(ui.actionFillet, SIGNAL(triggered()), this, SLOT(RandomExport()));
    connect(ui.actionChamfer, SIGNAL(triggered()), this, SLOT(MakePoint()));
    connect(ui.actionExtrude, SIGNAL(triggered()), this, SLOT(MakeCurve()));
    connect(ui.actionRevolve, SIGNAL(triggered()), this, SLOT(MakeSurface()));
    connect(ui.actionLoft, SIGNAL(triggered()), this, SLOT(ExportFile()));
    connect(ui.actionExtendCurveToPoint, SIGNAL(triggered()), this, SLOT(GenerateIsoCurves()));

    // Help
    connect(ui.actionAbout, SIGNAL(triggered()), this, SLOT(about()));
}

void occQt::createMenus(void)
{
}

void occQt::createToolBars(void)
{
    QToolBar* aToolBar;
    //aToolBar = addToolBar(tr("&Navigate"));
    //aToolBar->addAction(ui.actionZoom);
    //aToolBar->addAction(ui.actionPan);
    //aToolBar->addAction(ui.actionRotate);

    //aToolBar = addToolBar(tr("&View"));
    //aToolBar->addAction(ui.actionReset);
    //aToolBar->addAction(ui.actionFitAll);

    aToolBar = addToolBar(tr("&Intersection"));
    aToolBar->addAction(ui.actionBox);
    aToolBar->addAction(ui.actionCone);
    aToolBar->addAction(ui.actionSphere);
    aToolBar->addAction(ui.actionCylinder);
    aToolBar->addAction(ui.actionTorus);

    aToolBar = addToolBar(tr("&Utils"));
    aToolBar->addAction(ui.actionFillet);
    aToolBar->addAction(ui.actionChamfer);
    aToolBar->addAction(ui.actionExtrude);
    aToolBar->addAction(ui.actionRevolve);
    aToolBar->addAction(ui.actionLoft);
    aToolBar->addAction(ui.actionExtendCurveToPoint);
    aToolBar->addSeparator();
    aToolBar->addSeparator();

    //aToolBar = addToolBar(tr("About"));
    //aToolBar->addAction(ui.actionAbout);
}
Handle(Geom_BSplineCurve) IterateApproximate
(std::vector<Standard_Real>& InsertKnots,
    const std::vector<gp_Pnt>& Pnts,
    std::vector<Standard_Real>& PntsParams,
    std::vector<Standard_Real>& InitKnots,
    Standard_Integer degree,
    Standard_Integer MaxIterNum,
    Standard_Real toler);

std::vector<Standard_Real> ComputeUniformParam(Standard_Integer numSamples, Standard_Real left, Standard_Real right);
std::vector<Standard_Real> KnotGernerationByParams(const std::vector<Standard_Real>& params, Standard_Integer n, Standard_Integer p);
void createBSplineCurves()
{
    TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(MathTool::CreateCircleApproximation(10));
    BRepTools::Write(edge, "curve.brep");

    // 圆的参数
    Standard_Real radius = 10.0;    // 半径
    Standard_Real startAngle = 0.0; // 起始角度（0度）
    Standard_Real endAngle = 3 * M_PI / 2; // 四分之三圆 (270度)

    // 选择B-spline的控制点
    Standard_Integer numPoints = 10;  // 控制点的数量
    std::vector<gp_Pnt> intersectionPoints;
    // 随机数生成器（用于扰动）
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<Standard_Real> perturbation(-100, 100); // 设定扰动范围
    // 计算控制点

    for (Standard_Integer i = 0; i < numPoints; ++i)
    {
        Standard_Real angle = startAngle + i * (endAngle - startAngle) / (numPoints - 1);
        Standard_Real x;
        Standard_Real y;
        //if (i == numPoints / 2  - 1|| i == numPoints / 2 + 1)
        //{
        //    x = radius * cos(angle) * 3;
        //    y = radius * sin(angle) * -3;

        //}
        //else
        //{
        x = radius * cos(angle);
        y = radius * sin(angle);
        //}

        intersectionPoints.push_back(gp_Pnt(x, y, 0.0));
    }
    Handle(TColgp_HArray1OfPnt) points = new TColgp_HArray1OfPnt(1, intersectionPoints.size());
    for (Standard_Integer j = 0; j < intersectionPoints.size(); j++)
    {
        points->SetValue(j + 1, intersectionPoints[j]);
    }
    GeomAPI_Interpolate interpolate(points, Standard_False, 0.1);
    // 执行插值计算
    interpolate.Perform();
    // 检查是否成功完成插值
    //if (interpolate.IsDone())
    //{
    //    // 获取插值后的曲线对象
    //    Handle(Geom_BSplineCurve) aQuarterCircle = interpolate.Curve();
    //    TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(aQuarterCircle);
    //    BRepTools::Write(edge, "curve.brep");
    //}
    // 定义一系列点来创建折线段
    TColgp_Array1OfPnt LinePoints(1, 4); // 创建一个包含5个点的数组
    LinePoints.SetValue(1, gp_Pnt(0, 0, 0));
    LinePoints.SetValue(2, gp_Pnt(1, 2, 0));
    LinePoints.SetValue(3, gp_Pnt(2, 1, 0));
    LinePoints.SetValue(4, gp_Pnt(3, 4, 0));



    // 使用GeomAPI_PointsToBSpline将折线段转换为BSpline曲线
    GeomAPI_PointsToBSpline converter(LinePoints, 3, 9, GeomAbs_G2, 0.001);
    // 获取生成的BSpline曲线
    Handle(Geom_BSplineCurve) aLine = converter.Curve();
    edge = BRepBuilderAPI_MakeEdge(aLine);
    BRepTools::Write(edge, "aLine.brep");
}
std::vector<TopoDS_Edge> GetBSplineCurveD1(Handle(Geom_BSplineCurve)& theCurve, Standard_Integer theDiscreteNum)
{
    Standard_Real aFirstParameter = theCurve->FirstParameter();
    Standard_Real aLastParameter = theCurve->LastParameter();
    Standard_Real aStep = (aLastParameter - aFirstParameter) / theDiscreteNum;
    std::vector<TopoDS_Edge> aEdgeArray;
    for (Standard_Real t = aFirstParameter; t <= aLastParameter; t += aStep)
    {
        gp_Pnt aPnt;
        gp_Vec aVec;
        theCurve->D1(t, aPnt, aVec);
        aEdgeArray.push_back(BRepBuilderAPI_MakeEdge(aPnt, aPnt.Translated(aVec)));
    }
    return aEdgeArray;
}
std::vector<TopoDS_Edge> GetBSplineCurveD2(Handle(Geom_BSplineCurve)& theCurve, Standard_Integer theDiscreteNum)
{
    Standard_Real aFirstParameter = theCurve->FirstParameter();
    Standard_Real aLastParameter = theCurve->LastParameter();
    Standard_Real aStep = (aLastParameter - aFirstParameter) / theDiscreteNum;
    std::vector<TopoDS_Edge> aEdgeArray;
    for (Standard_Real t = aFirstParameter; t <= aLastParameter; t += aStep)
    {
        gp_Pnt aPnt;
        gp_Vec aVec1;
        gp_Vec aVec2;
        theCurve->D2(t, aPnt, aVec1, aVec2);
        aEdgeArray.push_back(BRepBuilderAPI_MakeEdge(aPnt, aPnt.Translated(aVec2)));
    }
    return aEdgeArray;
}
std::vector<TopoDS_Edge> GetBSplineCurveD3(Handle(Geom_BSplineCurve)& theCurve, Standard_Integer theDiscreteNum)
{
    Standard_Real aFirstParameter = theCurve->FirstParameter();
    Standard_Real aLastParameter = theCurve->LastParameter();
    Standard_Real aStep = (aLastParameter - aFirstParameter) / theDiscreteNum;
    std::vector<TopoDS_Edge> aEdgeArray;
    for (Standard_Real t = aFirstParameter; t <= aLastParameter; t += aStep)
    {
        gp_Pnt aPnt;
        gp_Vec aVec1;
        gp_Vec aVec2;
        gp_Vec aVec3;
        theCurve->D3(t, aPnt, aVec1, aVec2, aVec3);
        aEdgeArray.push_back(BRepBuilderAPI_MakeEdge(aPnt, aPnt.Translated(aVec3)));
    }
    return aEdgeArray;
}
//std::string format_as(const Handle(Geom_BSplineCurve) theCurve);
void UniformCurve(Handle(Geom_BSplineCurve)& curve);
std::vector<gp_Pnt> OccArrayConvertoVector(TColgp_Array1OfPnt theOccArray)
{
    std::vector<gp_Pnt> PntArray;
    for (int i = 1; i <= theOccArray.Size(); i++)
    {
        PntArray.push_back(theOccArray.Value(i));
    }
    return PntArray;
}
std::vector<gp_Pnt> OccArrayConvertoVector(Handle(Geom_BSplineCurve)& theCurve)
{
    TColgp_Array1OfPnt Poles = theCurve->Poles();
    std::vector<gp_Pnt> Result;
    for (int i = 1; i <= Poles.Size(); i++)
    {
        Result.push_back(Poles.Value(i));
    }
    return Result;
}
Standard_Real GetControlPointAverageDistance(std::vector<gp_Pnt> thePntArray)
{
    if (thePntArray.size() == 0) return 0;
    Standard_Real DistanceSum = 0;
    gp_Pnt LastPnt = thePntArray[0];
    for (int i = 1; i < thePntArray.size(); i++)
    {
        gp_Pnt CurrentPnt = thePntArray[i];
        DistanceSum += CurrentPnt.Distance(LastPnt);
    }
    return DistanceSum / thePntArray.size();
}

// 将 Standard_Real 转为保留两位小数的字符串
std::string format_as(Standard_Real value, int precision = 2) 
{
    std::ostringstream out;
    out << std::fixed << std::setprecision(precision) << value;
    return out.str();
}

// 读取 .dat 文件，返回 gp_Pnt 向量
std::vector<gp_Pnt> ReadDatToPoints(const std::string& filePath) 
{
    std::vector<gp_Pnt> points;
    std::ifstream infile(filePath);

    if (!infile.is_open())
    {
        return points;
    }

    std::string line;
    while (std::getline(infile, line)) 
    {
        std::istringstream iss(line);
        std::string x_str, y_str, z_str;

        if (std::getline(iss, x_str, ',') &&
            std::getline(iss, y_str, ',') &&
            std::getline(iss, z_str, ',')) 
        {

            try 
            {
                Standard_Real x = std::stod(x_str);
                Standard_Real y = std::stod(y_str);
                Standard_Real z = std::stod(z_str);
                points.emplace_back(x, y, z);
            }
            catch (const std::exception& e) 
            {
            }
        }
    }

    return points;
}

std::vector<std::string> GetFilesWithExtensions(
    const std::string& folderPath,
    const std::vector<std::string>& extensions,
    bool recursive = false)
{
    std::vector<std::string> result;

    try
    {
        if (recursive)
        {
            for (const auto& entry : std::filesystem::recursive_directory_iterator(folderPath))
            {
                if (!entry.is_regular_file()) continue;

                std::string ext = entry.path().extension().string();
                std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

                for (const auto& rawExt : extensions) 
                {
                    std::string normExt = rawExt;
                    if (!normExt.empty() && normExt[0] != '.')
                    {
                        normExt = "." + normExt;
                    }
                    std::transform(normExt.begin(), normExt.end(), normExt.begin(), ::tolower);

                    if (ext == normExt)
                    {
                        result.push_back(entry.path().string());
                        break;
                    }
                }
            }
        }
        else
        {
            for (const auto& entry : std::filesystem::directory_iterator(folderPath)) 
            {
                if (!entry.is_regular_file()) continue;

                std::string ext = entry.path().extension().string();
                std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

                for (const auto& rawExt : extensions) 
                {
                    std::string normExt = rawExt;
                    if (!normExt.empty() && normExt[0] != '.') 
                    {
                        normExt = "." + normExt;
                    }
                    std::transform(normExt.begin(), normExt.end(), normExt.begin(), ::tolower);

                    if (ext == normExt) 
                    {
                        result.push_back(entry.path().string());
                        break;
                    }
                }
            }
        }
    }
    catch (const std::exception& e) 
    {
        std::cerr << "遍历文件夹时出错: " << e.what() << std::endl;
    }

    return result;
}

Standard_Real AveragePointsDistance(const std::vector<gp_Pnt>& points, size_t K = 10) 
{
    if (points.size() < 2) return 0.0;

    gp_Pnt origin(0, 0, 0);

    std::vector<gp_Pnt> sortedPoints = points;
    std::sort(sortedPoints.begin(), sortedPoints.end(), [&origin](const gp_Pnt& a, const gp_Pnt& b) {
        return a.Distance(origin) < b.Distance(origin);
        });

    Standard_Real totalDist = 0.0;
    size_t count = 0;

    size_t n = sortedPoints.size();

    for (size_t i = 0; i < n; ++i) {
        size_t maxJ = std::min(i + K, n - 1);
        for (size_t j = i + 1; j <= maxJ; ++j) {
            totalDist += sortedPoints[i].Distance(sortedPoints[j]);
            ++count;
        }
    }

    return count == 0 ? 0.0 : totalDist / count;
}
void occQt::about()
{
    // 弹出文件夹选择对话框
    QString folderPath = QFileDialog::getExistingDirectory(nullptr, "选择文件夹", QDir::homePath());

    if (folderPath.isEmpty())
    {
        return;
    }

    QDir folderDir(folderPath);


    //std::vector<std::string> filePaths = GetFilesWithExtensions(folderPath.toStdString(), {".dat", ".brep"});
    //std::vector<std::string> filePaths = GetFilesWithExtensions(folderPath.toStdString(), {".brep"});
    std::vector<std::string> filePaths = GetFilesWithExtensions(folderPath.toStdString(), {".dat"},false);
    std::cout << filePaths.size() << std::endl;
    std::string doneDir = folderPath.toStdString() + "/done";
    std::string errorDir = folderPath.toStdString() + "/error";
    std::filesystem::create_directories(doneDir);
    std::filesystem::create_directories(errorDir);

    for (auto file : filePaths)
    {
        bool success = false;
        auto dot = file.find_last_of('.');
        auto fileName = file.substr(0, dot);
        auto extension = file.substr(dot + 1);
        std::string outputFolder = fileName + "_Result";
        std::string fileOnlyName = std::filesystem::path(file).filename().string();
        std::vector<gp_Pnt> points;
        Handle(Geom_BSplineCurve) curve;

        try
        {
            // 清空已有输出文件夹内容
            if (std::filesystem::exists(outputFolder))
            {
                for (const auto& entry : std::filesystem::directory_iterator(outputFolder))
                {
                    std::filesystem::remove_all(entry.path());
                }
            }
            else
            {
                std::filesystem::create_directory(outputFolder);
            }

            Standard_Real PointAverageDistance;
            CurveFair::ExportFilePath = outputFolder;

            if (extension == "dat")
            {
                CurveFair::ExportFilePath = outputFolder + "/FitPointResult";
                std::filesystem::create_directory(CurveFair::ExportFilePath);

                points = ReadDatToPoints(file);
                std::cout << "Points size:" << points.size() << std::endl;
                PointAverageDistance = AveragePointsDistance(points);
                Visualize(points, Quantity_NOC_BLUE);

                Handle(Geom_BSplineCurve) originalCurve = CurveFair::SampleAndFitBSpline(points, FITTOLERANCE);
                std::vector<Handle(Geom_BSplineCurve)> resultSegments = { originalCurve };
                std::vector<std::vector<gp_Pnt>> fitPointsSegments = { points };
                SurfaceModelingTool::ExportPoints(points, CurveFair::ExportFilePath + "/FitPoints.step");
                SurfaceModelingTool::ExportBSplineCurves({ originalCurve }, CurveFair::ExportFilePath + "/OriginalCurve.step");

                //Standard_Real totalLength = GCPnts_AbscissaPoint::Length(GeomAdaptor_Curve(originalCurve));
                //Standard_Real segmentMaxLength = 30000;
                //Standard_Integer nSegments = std::ceil(totalLength / segmentMaxLength);


                //fitPointsSegments.resize(nSegments);
                //for (const auto& pt : points)
                //{
                //    // 最近点投影求 t 参数
                //    GeomAPI_ProjectPointOnCurve projector(pt, originalCurve);
                //    if (!projector.NbPoints()) continue;

                //    Standard_Real t = projector.LowerDistanceParameter();

                //    // 弧长计算（从曲线起点到投影点）
                //    Standard_Real s = GCPnts_AbscissaPoint::Length(GeomAdaptor_Curve(originalCurve), originalCurve->FirstParameter(), t);

                //    int idx = std::min(int(s / segmentMaxLength), nSegments - 1);
                //    fitPointsSegments[idx].push_back(pt);
                //}

                //Standard_Real curveFirstParam = originalCurve->FirstParameter();
                //Standard_Real curveLastParam = originalCurve->LastParameter();
                //Standard_Real tolerance = 1e-6;  // 你可以根据需求调整

                //std::vector<Handle(Geom_BSplineCurve)> resultSegments;
                //for (int i = 0; i < nSegments; ++i)
                //{
                //    Standard_Real s_start = i * segmentMaxLength;
                //    Standard_Real s_end = std::min((i + 1) * segmentMaxLength, totalLength);

                //    Standard_Real t_start = (i == 0) ? curveFirstParam :
                //        MathTool::GetParameterAtLength(originalCurve, curveFirstParam, s_start, tolerance);
                //    Standard_Real t_end = (i == nSegments - 1) ? curveLastParam :
                //        MathTool::GetParameterAtLength(originalCurve, curveFirstParam, s_end, tolerance);
                //    // 用 Geom_TrimmedCurve 截取曲线段
                //    Handle(Geom_TrimmedCurve) segment = new Geom_TrimmedCurve(originalCurve, t_start, t_end);


                //    // CurveFair 的构造函数需要 Geom_BSplineCurve，因此要转换
                //    Handle(Geom_BSplineCurve) segmentBSpline = GeomConvert::CurveToBSplineCurve(segment);;
                //    UniformCurve(segmentBSpline);
                //    resultSegments.push_back(segmentBSpline);

                //}
                for (Standard_Real hausdorffDistance = 50; hausdorffDistance >= 50; hausdorffDistance *= 0.5)
                {
                    Standard_Real hausdorffDistanceResult = INT_MIN;
                    std::vector<Handle(Geom_BSplineCurve)> results;
                    for (int i = 0; i < resultSegments.size(); i++)
                    {
                        std::cout << "处理前节点数量" << resultSegments[i]->NbKnots() << std::endl;
                        std::cout << "处理前控制点数量" << resultSegments[i]->NbPoles() << std::endl;
                        CurveFair aCurveFairConstructor(resultSegments[i], fitPointsSegments[i], hausdorffDistance);
                        results.push_back(aCurveFairConstructor.GetResult());
                        Visualize(aCurveFairConstructor.GetOriginalCurve(), Quantity_NOC_GOLD);
                        std::cout << "处理后节点数量" <<aCurveFairConstructor.GetOriginalCurve()->NbKnots() << std::endl;
                        std::cout << "处理后控制点数量" <<aCurveFairConstructor.GetOriginalCurve()->NbPoles() << std::endl;
                        Visualize(aCurveFairConstructor.GetResult(), Quantity_NOC_RED);
                        Visualize(OccArrayConvertoVector(aCurveFairConstructor.GetResult()->Poles()), Quantity_NOC_RED);
                        Standard_Real hsResult = aCurveFairConstructor.GetHausdorffDistanceResult();
                        hausdorffDistanceResult = std::max(hsResult, hausdorffDistanceResult);

                    }

                    GeomConvert_CompCurveToBSplineCurve joiner(results[0]);
                    Standard_Boolean isJoint = Standard_True;
                    for (Standard_Integer i = 1; i < results.size(); ++i)
                    {
                        if (!joiner.Add(results[i], Precision::Confusion()))
                        {
                            std::cout << "Failed to join curve at index " << i << std::endl;
                            isJoint = Standard_False;
                            break;
                        }
                    }

                    if (!isJoint)
                    {
                        throw("ERROR::FAILED TO JOIN CURVES!!!");
                    }
                    SurfaceModelingTool::ExportBSplineCurves({ curve = joiner.BSplineCurve() },
                        CurveFair::ExportFilePath + "/Result_HausdorffDistance_" + format_as(hausdorffDistanceResult) + ".step");
                }

                myOccView->getContext()->EraseAll(Standard_True);
            }
            else if (extension == "brep")
            {
                CurveFair::ExportFilePath = outputFolder + "/CurveResult";
                std::filesystem::create_directory(outputFolder);
                std::vector<Handle(Geom_BSplineCurve)> aCurveArray;
                SurfaceModelingTool::LoadBSplineCurves(file, aCurveArray);

                GeomConvert_CompCurveToBSplineCurve joiner(aCurveArray[0]);
                Standard_Boolean isJoint = Standard_True;
                for (Standard_Integer i = 1; i < aCurveArray.size(); ++i)
                {
                    if (!joiner.Add(aCurveArray[i], Precision::Confusion()))
                    {
                        std::cout << "Failed to join curve at index " << i << std::endl;
                        isJoint = Standard_False;
                        break;
                    }
                }

                if (!isJoint)
                {
                    throw("ERROR::FAILED TO JOIN CURVES!!!");
                }

                curve = joiner.BSplineCurve();
                PointAverageDistance = GetControlPointAverageDistance(OccArrayConvertoVector(curve));
                for (Standard_Real hausdorffDistance = 1000; hausdorffDistance >= 30; hausdorffDistance *= 0.8)
                {
                    CurveFair aCurveFairConstructor(curve, points, hausdorffDistance);
                    gp_Vec tangent1(aCurveFairConstructor.GetOriginalCurve()->Pole(1), aCurveFairConstructor.GetOriginalCurve()->Pole(2));
                    gp_Vec tangent2(aCurveFairConstructor.GetOriginalCurve()->Pole(aCurveFairConstructor.GetOriginalCurve()->NbPoles() - 1), aCurveFairConstructor.GetOriginalCurve()->Pole(aCurveFairConstructor.GetOriginalCurve()->NbPoles()));

                    SurfaceModelingTool::ExportBSplineCurves({ curve }, CurveFair::ExportFilePath + "/OriginalCurve.step");

                    gp_Vec newTangent1(aCurveFairConstructor.GetResult()->Pole(1), aCurveFairConstructor.GetResult()->Pole(2));
                    gp_Vec newTangent2(aCurveFairConstructor.GetResult()->Pole(aCurveFairConstructor.GetResult()->NbPoles() - 1), aCurveFairConstructor.GetResult()->Pole(aCurveFairConstructor.GetResult()->NbPoles()));

                    Standard_Real angle1 = tangent1.Angle(newTangent1);
                    Standard_Real angle2 = tangent2.Angle(newTangent2);
                    std::cout << "ERROR::Tangent vectors are not equal." << "Angle1 : " << angle1 << " Angle2 : " << angle2 << std::endl;

                    Visualize(aCurveFairConstructor.GetOriginalCurve(), Quantity_NOC_GOLD);
                    Visualize(aCurveFairConstructor.GetResult(), Quantity_NOC_RED); 
                    Visualize(OccArrayConvertoVector(aCurveFairConstructor.GetResult()->Poles()), Quantity_NOC_RED);
                    Standard_Real hausdorffDistanceResult = aCurveFairConstructor.GetHausdorffDistanceResult();
                    SurfaceModelingTool::ExportBSplineCurves({ aCurveFairConstructor.GetResult() }, CurveFair::ExportFilePath + "/Result_HausdorffDistance_" + format_as(hausdorffDistanceResult) + ".step");
                }
                myOccView->getContext()->EraseAll(Standard_True);
            }

            success = true;
        }
        catch (const std::exception& e)
        {
            std::cerr << "处理文件出错: " << e.what() << std::endl;
        }

        // 成功或失败后移动原始文件
        try
        {
            std::string destination = success ? (doneDir + "/" + fileOnlyName) : (errorDir + "/" + fileOnlyName);
            std::filesystem::rename(file, destination);
        }
        catch (const std::exception& e)
        {
            std::cerr << "移动文件失败: " << e.what() << std::endl;
        }
    }
}


void occQt::importFile()
{
    myOccView->getContext()->EraseAll(Standard_True);
    TopoDS_Shape aTopoBox = BRepPrimAPI_MakeBox(3.0, 4.0, 5.0).Shape();
    TopoDS_Shape shape;
    STEPControl_Reader reader;
    QString filename = QFileDialog::getOpenFileName(this, "选择文件", QDir::homePath(), "All Files (*.*);;Text Files (*.txt);;Image Files (*.png *.jpg)");
    if (reader.ReadFile(filename.toStdString().c_str()) == IFSelect_RetDone) {
        reader.TransferRoots();
        shape = reader.OneShape();
    }
    TopExp_Explorer exp;
    Standard_Integer cnt = 0;
    for (exp.Init(shape, TopAbs_FACE); exp.More(); exp.Next())
    {
        if (!cnt) {
            shape1 = exp.Current();
        }
        else
        {
            shape2 = exp.Current();
        }
        cnt++;
    }
    //Handle(AIS_Shape) anAisBox = new AIS_Shape(aTopoBox);
    Handle(AIS_Shape) aShape1 = new AIS_Shape(shape1);
    Handle(AIS_Shape) aShape2 = new AIS_Shape(shape2);

    //anAisBox->SetColor(Quantity_NOC_AZURE);
    aShape1->SetColor(Quantity_NOC_BISQUE);
    aShape2->SetColor(Quantity_NOC_BISQUE);

    myOccView->getContext()->Display(aShape1, Standard_True);
    myOccView->getContext()->Display(aShape2, Standard_True);
    intersectionDone = Standard_False;
}

void occQt::Triangulation()
{
    std::vector<TopoDS_Shape> shapes;
    Standard_Integer choice;

    if (shape1.IsNull() && shape2.IsNull())
    {
        QString filename = QFileDialog::getOpenFileName(this, "choose files", QDir::homePath(), "All Files (*.*);;Text Files (*.txt);;Image Files (*.png *.jpg)");
        ModelImporter importer;
        shapes = importer.loadStepShape(filename.toStdString());
        shape1 = shapes[0];
        shape2 = shapes[1];
    }
    else
    {
        shapes.push_back(shape1);
        shapes.push_back(shape2);
    }

    // Display shape1 and shape2
    myOccView->getContext()->RemoveAll(Standard_True);
    Handle(AIS_Shape) aShape1 = new AIS_Shape(shape1);
    aShape1->SetColor(Quantity_NOC_BISQUE);
    myOccView->getContext()->Display(aShape1, Standard_True);
    Handle(AIS_Shape) aShape2 = new AIS_Shape(shape2);
    aShape2->SetColor(Quantity_NOC_BISQUE);
    myOccView->getContext()->Display(aShape2, Standard_True);

    Standard_Integer cnt = 0;
    // Store triangulated results and build R-trees
    for (auto shape : shapes)
    {
        faceList[cnt] = ModelTriangulation::ModelTriangulate(shape, triPoints[cnt]);
        cnt++;
    }

    choice = QMessageBox::question(this, "Confirmation", "Display Results? ", QMessageBox::Yes | QMessageBox::No);
    if (choice == QMessageBox::Yes)
    {
        myOccView->getContext()->RemoveAll(Standard_True);
        std::vector<std::future<void>> futures;
        std::mutex displayMutex;

        for (Standard_Integer j = 0; j <= 1; j++)
        {
            for (size_t i = 0; i < faceList[j].size(); ++i)
            {
                futures.push_back(std::async(std::launch::async, [&, j, i]()
                    {
                        Handle(AIS_Shape) triShape = new AIS_Shape(faceList[j][i]);
                        triShape->SetColor(Quantity_NOC_BISQUE);
                        // 互斥锁,析构时自动解锁
                        std::lock_guard<std::mutex> lock(displayMutex);
                        myOccView->getContext()->Display(triShape, Standard_True);
                    }));
            }
        }

        // 等待所有异步任务完成
        for (auto& f : futures) {
            f.get();
        }
    }
}


void occQt::TriangulationIntersection()
{
    if (faceList[0].size() && faceList[1].size())
    {
        auto start_time1 = std::chrono::high_resolution_clock::now();
        // 创建R树
        R_Tree RTree1, RTree2;
        RTree1.BuildRtree(faceList[0], triPoints[0]);
        RTree2.BuildRtree(faceList[1], triPoints[1]);
        // 获取结束时间点

        auto end_time1 = std::chrono::high_resolution_clock::now();
        // 计算程序执行时间，单位为秒
        std::chrono::duration<Standard_Real> execution_time1 = end_time1 - start_time1;

        // 将执行时间转换为毫秒
        Standard_Real execution_time_ms1 = execution_time1.count() * 1000.0;

        std::cout << "程序执行时间为：" << execution_time_ms1 << " 毫秒" << std::endl;
        // 初始化结果
        result.Nullify();
        BRep_Builder builder;
        builder.MakeCompound(result);

        std::vector<std::future<void>> futures;
        std::mutex sectionMutex; // 用于保护BRepAlgoAPI_Section操作

        // 定义一个递归函数来处理每个节点
        std::function<void(R_TreeNode*)> processNode;
        processNode = [&](R_TreeNode* currentNode) {
            // 在第二个 R 树中找到可能相交的节点
            std::vector<R_TreeNode*> potentialIntersectNodes = RTree2.Search(currentNode->Box, RTree2.Root);

            for (R_TreeNode* node2 : potentialIntersectNodes)
            {
                // 将每个求交操作作为一个异步任务
                futures.push_back(std::async(std::launch::async, [&, currentNode, node2]()
                    {
                        std::unique_lock<std::mutex> lock(sectionMutex);
                        BRepAlgoAPI_Section section(currentNode->Face, node2->Face, Standard_True);
                        section.ComputePCurveOn1(Standard_True);
                        section.Approximation(Standard_True);
                        section.Build();
                        lock.unlock(); // 释放锁

                        if (section.IsDone())
                        {
                            TopoDS_Shape intersectionLine = section.Shape();
                            // 在这里不需要锁，因为builder.Add被假定为线程安全或对结果的影响不会造成冲突
                            builder.Add(result, intersectionLine);
                        }
                    }));
            }

            // 将当前节点的子节点加入递归处理
            for (R_TreeNode* child : currentNode->Childs)
            {
                processNode(child);
            }
            };

        // 开始从根节点处理
        processNode(RTree1.Root);

        // 等待所有异步任务完成
        for (auto& future : futures) {
            future.get();
        }

        Handle(AIS_Shape) aIntersectionLine = new AIS_Shape(result);
        aIntersectionLine->SetColor(Quantity_NOC_RED);
        myOccView->getContext()->Display(aIntersectionLine, Standard_True);
        faceList[0].clear();
        faceList[1].clear();
        intersectionDone = Standard_True;
    }
    else
    {
        Triangulation();
        TriangulationIntersection();
    }

}

void occQt::ClearDisplay()
{
    myOccView->getContext()->RemoveAll(Standard_True);
    shape1 = shape2 = TopoDS_Shape();
    curveShape = TopoDS_Shape();
    pointShape = TopoDS_Shape();
    faceList[0].clear();
    faceList[1].clear();
    result.Nullify();
    boundedCurve.Nullify();

}

void occQt::PrintInfo()
{
    // 在文本框中输出顶点信息
    QString outputText;
    TopExp_Explorer exp;
    Standard_Integer cnt = 1;
    for (exp.Init(result, TopAbs_EDGE); exp.More(); exp.Next()) {
        TopoDS_Edge edge = TopoDS::Edge(exp.Current());
        // 获取边的起始和结束顶点
        TopoDS_Vertex startVertex, endVertex;
        TopExp::Vertices(edge, startVertex, endVertex);

        // 获取顶点的坐标
        gp_Pnt startPoint = BRep_Tool::Pnt(startVertex);
        gp_Pnt endPoint = BRep_Tool::Pnt(endVertex);

        // 将顶点信息添加到输出文本
        outputText += QString("Edge%1: startV : (%2, %3, %4), endV : (%5, %6, %7)\n")
            .arg(cnt).arg(startPoint.X()).arg(startPoint.Y()).arg(startPoint.Z())
            .arg(endPoint.X()).arg(endPoint.Y()).arg(endPoint.Z());
        cnt++;
    }

    // 在对话框中显示文本框
    QDialog dialog(this);
    QVBoxLayout* layout = new QVBoxLayout(&dialog);


    QTextEdit* outputTextEdit = new QTextEdit(&dialog);
    outputTextEdit->setPlainText(outputText);
    // 设置文本框对象名称
    layout->addWidget(outputTextEdit);

    // 设置文本框扩展以填充空间
    outputTextEdit->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    dialog.setWindowTitle("Edge Info");

    // 显示对话框
    dialog.exec();

}

void occQt::RandomExport()
{

    // 创建一个文件对话框
    QFileDialog dialog;

    // 设置对话框模式为选择文件夹
    dialog.setFileMode(QFileDialog::Directory);
    QString folderPath;
    // 显示对话框
    if (dialog.exec()) {
        // 获取用户选择的文件夹路径
        folderPath = dialog.selectedFiles().first();
    }
    // 获取两个简单字符串
    bool ok;
    QString string1 = QInputDialog::getText(nullptr, "Input", "Enter String 1: (Basic Model Type: plane, cone, cylinder, sphere, torus)", QLineEdit::Normal, "", &ok);
    QString string2 = QInputDialog::getText(nullptr, "Input", "Enter String 2: (Basic Model Type: plane, cone, cylinder, sphere, torus)", QLineEdit::Normal, "", &ok);
    Standard_Integer value = QInputDialog::getInt(nullptr, "Input", "Enter an Integer:", 0, 0, 100, 1, &ok);


    // 使用用户输入调用RandomExport::randomRotateAndExport()
    RandomExport::randomRotateAndExport(folderPath.toStdString().c_str(), string1.toStdString().c_str(), string2.toStdString().c_str(), value);
}


void occQt::MakePoint()
{
    bool ok;

    Standard_Real x = QInputDialog::getDouble(nullptr, "Input", "Enter X for point:", 0, -1000, 1000, 8, &ok);
    Standard_Real y = QInputDialog::getDouble(nullptr, "Input", "Enter Y for point:", 0, -1000, 1000, 8, &ok);
    Standard_Real z = QInputDialog::getDouble(nullptr, "Input", "Enter Z for point:", 0, -1000, 1000, 8, &ok);

    TopoDS_Shape point = MakeShape::MakePoint(x, y, z);
    if (!point.IsNull())
    {
        Handle(AIS_Shape) aPoint = new AIS_Shape(point);
        aPoint->SetColor(Quantity_NOC_BISQUE);
        myOccView->getContext()->Display(aPoint, Standard_True);
    }
}


void occQt::MakeCurve()
{
    bool ok;
    QString curveType = QInputDialog::getText(nullptr, "Input", "Enter Curve Type: (circle, ellipse, parabola, hyperbola, line)", QLineEdit::Normal, "", &ok);
    if (ok)
        curveShape = createCurve(curveType);

    if (!curveShape.IsNull())
    {
        Handle(AIS_Shape) aCurve = new AIS_Shape(curveShape);
        aCurve->SetColor(Quantity_NOC_BISQUE);
        myOccView->getContext()->Display(aCurve, Standard_True);
    }
}

void occQt::MakeSurface()
{
    bool ok;
    TopoDS_Shape shape;
    QString surfaceType = QInputDialog::getText(nullptr, "Input", "Enter Surface Type: (plane, cylinder, sphere, cone, torus)", QLineEdit::Normal, "", &ok);
    if (ok)
    {
        shape = createSurface(surfaceType);
        if (changeCount % 2 == 0)
        {
            shape1 = shape;
        }
        else
        {
            shape2 = shape;
        }
        changeCount++;
    }

    if (!shape.IsNull())
    {
        Handle(AIS_Shape) aShape = new AIS_Shape(shape);
        aShape->SetColor(Quantity_NOC_BISQUE);
        myOccView->getContext()->Display(aShape, Standard_True);
    }
    intersectionDone = Standard_False;
}

void occQt::ExportFile()
{
    // 创建一个文件对话框
    QFileDialog dialog;

    // 设置对话框模式为选择文件夹
    dialog.setFileMode(QFileDialog::Directory);
    QString folderPath;
    // 显示对话框
    if (dialog.exec()) {
        // 获取用户选择的文件夹路径
        folderPath = dialog.selectedFiles().first();
    }

    // 指定文件名
    QString filePath = folderPath + "/exportModel.step";

    // 获取文件的基本信息
    QFileInfo fileInfo(filePath);

    // 如果文件已存在，则自动添加后缀
    Standard_Integer suffix = 1;
    while (fileInfo.exists()) {
        filePath = folderPath + QString("/exportModel_%1.step").arg(suffix);
        fileInfo.setFile(filePath);
        suffix++;
    }

    BRep_Builder builder;
    TopoDS_Compound exportResult;
    builder.MakeCompound(exportResult);
    if (intersectionDone)
    {
        builder.Add(exportResult, result);
    }
    builder.Add(exportResult, shape1);
    builder.Add(exportResult, shape2);

    ModelExporter exporter;
    exporter.exportModel(exportResult, filePath.toStdString());

    exportResult.Nullify();

}

// 对单条等参线进行等距采样并计算弧长
std::vector<std::pair<gp_Pnt, Standard_Real>> SampleIsoCurveWithArcLength(const Handle(Geom_BSplineCurve)& bsplineCurve, Standard_Integer numSamples)
{
    Handle(GeomAdaptor_Curve) curve = new GeomAdaptor_Curve(bsplineCurve);
    GCPnts_QuasiUniformAbscissa sampler(*curve, numSamples);
    std::vector<std::pair<gp_Pnt, Standard_Real>> sampledPointsWithArcLength;

    for (Standard_Integer i = 1; i <= sampler.NbPoints(); ++i)
    {
        Standard_Real param = sampler.Parameter(i);
        gp_Pnt point = curve->Value(param); // 获取点坐标

        // 计算当前点与第一个点之间的弧长
        Standard_Real arcLength = CPnts_AbscissaPoint::Length(*curve, sampler.Parameter(1), param);

        // 更新点和对应弧长存储
        sampledPointsWithArcLength.emplace_back(point, arcLength);
    }

    return sampledPointsWithArcLength;
}
// 输出采样点和对应弧长的函数
void PrintSampledPointsWithArcLength(const std::vector<std::pair<gp_Pnt, Standard_Real>>& sampledPoints) {
    for (const auto& pair : sampledPoints)
    {
        const gp_Pnt& point = pair.first;
        Standard_Real arcLength = pair.second;
        std::cout << "Point: (" << point.X() << ", " << point.Y() << ", " << point.Z() << "), Arc Length: " << arcLength << std::endl;
    }
}

Handle(Geom_BSplineSurface) GenerateCoonsSurface(
    Handle(Geom_BSplineCurve)& curve1, Handle(Geom_BSplineCurve)& curve2, Handle(Geom_BSplineCurve)& curve3, Handle(Geom_BSplineCurve)& curve4
) {

    // 创建 GeomFill_BSplineCurves 对象，使用 Coons 填充样式
    GeomFill_BSplineCurves fillCoonsStyle(
        curve1, // U=0
        curve2, // V=1
        curve3, // U=1
        curve4, // V=0
        GeomFill_CoonsStyle
    );

    // 获取生成的曲面
    Handle(Geom_Surface) surface = fillCoonsStyle.Surface();

    // 检查曲面是否生成成功
    if (surface.IsNull()) {
        throw std::runtime_error("使用 CoonsStyle 创建曲面失败！");
    }

    std::cout << "使用 CoonsStyle 成功创建曲面！" << std::endl;

    // 尝试将 Geom_Surface 转换为 Geom_BSplineSurface
    Handle(Geom_BSplineSurface) bsplineSurface = Handle(Geom_BSplineSurface)::DownCast(surface);
    if (bsplineSurface.IsNull())
    {
        throw std::runtime_error("生成的曲面不是 B-Spline 曲面！");
    }

    return bsplineSurface;
}

//Uniform Curve Knot to [0,1]
void UniformCurve(Handle(Geom_BSplineCurve)& curve);

void occQt::GenerateIsoCurves(void)
{

    for (Standard_Integer i = 35; i <= 53; i++)
    {
        if (i == 23 || (i >= 42 && i <= 45) || i == 40) continue;;
        occQt::isVisualize = Standard_True;
        myOccView->getContext()->RemoveAll(Standard_True);
        // 读入边界线
        std::vector<Handle(Geom_BSplineCurve)> tempArray;
        tempArray.clear();
        std::string filename = "E:/Models/";
        filename += std::to_string(i);
        filename += "_";
        std::string boundaryPath = filename + "boundary.brep";
        SurfaceModelingTool::LoadBSplineCurves(boundaryPath.c_str(), tempArray);
        SurfaceModelingTool::ApproximateBoundaryCurves(tempArray);

        Visualize(tempArray);
        for (Standard_Integer j = 0; j < tempArray.size(); j++)
        {
            Standard_Real x = (tempArray[j]->StartPoint().X() + tempArray[j]->EndPoint().X()) / 2;
            Standard_Real y = (tempArray[j]->StartPoint().Y() + tempArray[j]->EndPoint().Y()) / 2;
            Standard_Real z = (tempArray[j]->StartPoint().Z() + tempArray[j]->EndPoint().Z()) / 2;
            gp_Pnt pnt(x, y, z);
            Standard_Real value = Standard_Real(j + 1);
            Visualize(std::make_pair(pnt, value));
        }

        std::vector<Handle(Geom_BSplineCurve)> anInternalBSplineCurves;
        std::string internalPath = filename + "internal.brep";

        SurfaceModelingTool::LoadBSplineCurves(internalPath, anInternalBSplineCurves);
        Visualize(anInternalBSplineCurves);
        MathTool::TrimInternalCurves(anInternalBSplineCurves, tempArray);
        Visualize(anInternalBSplineCurves);
        if (tempArray.size() == 3)
        {
            // 初始化一个向量用于存储每条曲线的交点计数以及对应的样条曲线
            std::vector<std::pair<Standard_Integer, Handle(Geom_BSplineCurve)>> anInterCount =
            {
                {0, tempArray[0]},
                {0, tempArray[1]},
                {0, tempArray[2]}
            };

            // 遍历内部样条曲线与tempArray中的曲线，统计每条曲线的交点计数
            for (auto internalCurve : anInternalBSplineCurves)
            {
                for (Standard_Integer m = 0; m < tempArray.size(); m++)
                {
                    if (MathTool::ComputeAngleBetweenCurves(tempArray[m], internalCurve) < 5.0)
                    {
                        anInterCount[m].first++;
                    }
                }
            }

            // 按交点计数从大到小对曲线进行排序
            std::sort(anInterCount.begin(), anInterCount.end(),
                [](const std::pair<Standard_Integer, Handle(Geom_BSplineCurve)>& curve1,
                    const std::pair<Standard_Integer, Handle(Geom_BSplineCurve)>& curve2)
                {
                    return curve1.first < curve2.first;
                });

            if (anInterCount[2].first >= 100)
            {
                // 清空并重新调整tempArray，将交点最多的两条边作为u方向的边
                tempArray.clear();
                tempArray.resize(4);
                tempArray[0] = anInterCount[0].second;
                tempArray[2] = anInterCount[1].second;
                tempArray[1] = anInterCount[2].second;

                // 将交点最多的两条边的交点作为退化边
                gp_Pnt aDegeneratePoint(0, 0, 0);
                if (tempArray[0]->StartPoint().Distance(tempArray[2]->StartPoint()) > 10
                    && tempArray[0]->StartPoint().Distance(tempArray[2]->EndPoint()) > 10)
                {
                    aDegeneratePoint = tempArray[0]->EndPoint();
                }
                else if (tempArray[0]->EndPoint().Distance(tempArray[2]->StartPoint()) > 10
                    && tempArray[0]->EndPoint().Distance(tempArray[2]->EndPoint()) > 10)
                {
                    aDegeneratePoint = tempArray[0]->StartPoint();
                }

                // 构建退化边
                TColgp_Array1OfPnt poles(1, 2);
                poles.SetValue(1, aDegeneratePoint);
                poles.SetValue(2, aDegeneratePoint);

                TColStd_Array1OfReal knots(1, 2);
                knots.SetValue(1, 0.0);
                knots.SetValue(2, 1.0);

                TColStd_Array1OfInteger multiplicities(1, 2);
                multiplicities.SetValue(1, 2);
                multiplicities.SetValue(2, 2);

                tempArray[3] = new Geom_BSplineCurve(poles, knots, multiplicities, 1);

                // 调整次序，要求首尾相接
                Standard_Real tol = tempArray[0]->EndPoint().Distance(tempArray[0]->StartPoint()) / 1000;
                if (tempArray[0]->StartPoint().Distance(tempArray[1]->StartPoint()) < tol)
                {
                    tempArray[0]->Reverse();
                }
                if (tempArray[2]->EndPoint().Distance(tempArray[1]->EndPoint()) < tol)
                {
                    tempArray[2]->Reverse();
                }
            }
            else
            {
                gp_Pnt pnt1 = tempArray[0]->StartPoint(), pnt2 = tempArray[0]->EndPoint(), pnt3 = tempArray[1]->StartPoint();
                gp_Pnt pnt4 = tempArray[1]->EndPoint(), pnt5 = tempArray[2]->StartPoint(), pnt6 = tempArray[2]->EndPoint();

                Standard_Real tol = tempArray[0]->EndPoint().Distance(tempArray[0]->StartPoint()) / 1000;
                if (tempArray[1]->EndPoint().Distance(tempArray[0]->EndPoint()) < tol)
                {
                    tempArray[1]->Reverse();
                }
                else if (tempArray[2]->StartPoint().Distance(tempArray[0]->EndPoint()) < tol)
                {
                    std::swap(tempArray[1], tempArray[2]);
                }
                else if (tempArray[2]->EndPoint().Distance(tempArray[0]->EndPoint()) < tol)
                {
                    std::swap(tempArray[1], tempArray[2]);
                    tempArray[1]->Reverse();
                }

                if (tempArray[2]->EndPoint().Distance(tempArray[1]->EndPoint()) < tol)
                {
                    tempArray[2]->Reverse();
                }

                pnt1 = tempArray[0]->StartPoint(); pnt2 = tempArray[0]->EndPoint(); pnt3 = tempArray[1]->StartPoint();
                pnt4 = tempArray[1]->EndPoint(); pnt5 = tempArray[2]->StartPoint(); pnt6 = tempArray[2]->EndPoint();

                // 三边情况，创建退化边构成四边
                std::vector<gp_Pnt> boundaryPoints = { tempArray[0]->StartPoint(), tempArray[1]->StartPoint(), tempArray[2]->StartPoint() };

                // 定义边
                gp_Vec line_01(boundaryPoints[1].XYZ() - boundaryPoints[0].XYZ());
                gp_Vec line_12(boundaryPoints[2].XYZ() - boundaryPoints[1].XYZ());
                gp_Vec line_20(boundaryPoints[0].XYZ() - boundaryPoints[2].XYZ());

                auto calculateAngle = [](const gp_Vec& v1, const gp_Vec& v2)
                    {
                        Standard_Real dotProduct = v1.Dot(v2);
                        Standard_Real magnitudes = v1.Magnitude() * v2.Magnitude();
                        return std::acos(dotProduct / magnitudes);  // 返回角度
                    };

                // 计算三个角的夹角
                Standard_Real angleAtPoint0 = calculateAngle(-line_20, line_01);  // 点0的夹角
                Standard_Real angleAtPoint1 = calculateAngle(-line_01, line_12);  // 点1的夹角
                Standard_Real angleAtPoint2 = calculateAngle(-line_12, line_20);  // 点2的夹角

                // 找出最大角度
                Standard_Real maxAngle = std::max({ angleAtPoint0, angleAtPoint1, angleAtPoint2 });
                Standard_Integer maxAngleIndex = 0;
                if (maxAngle == angleAtPoint1) maxAngleIndex = 1;
                else if (maxAngle == angleAtPoint2) maxAngleIndex = 2;

                //// 构造退化边
                //auto CreateDegenerateEdge = [](const gp_Pnt& point)
                //    {
                //        TColgp_Array1OfPnt poles(1, 2);
                //        poles.SetValue(1, point);
                //        poles.SetValue(2, point);

                //        TColStd_Array1OfReal knots(1, 2);
                //        knots.SetValue(1, 0.0);
                //        knots.SetValue(2, 1.0);

                //        TColStd_Array1OfInteger multiplicities(1, 2);
                //        multiplicities.SetValue(1, 2);
                //        multiplicities.SetValue(2, 2);

                //        return new Geom_BSplineCurve(poles, knots, multiplicities, 1);
                //    };

                // 构造退化边
                auto CreateDegenerateEdge = [](const gp_Pnt& p1, const gp_Pnt& p2)
                    {
                        TColgp_Array1OfPnt poles(1, 2);
                        poles.SetValue(1, p1);
                        poles.SetValue(2, p2);

                        TColStd_Array1OfReal knots(1, 2);
                        knots.SetValue(1, 0.0);
                        knots.SetValue(2, 1.0);

                        TColStd_Array1OfInteger multiplicities(1, 2);
                        multiplicities.SetValue(1, 2);
                        multiplicities.SetValue(2, 2);

                        return new Geom_BSplineCurve(poles, knots, multiplicities, 1);
                    };
                if (maxAngleIndex == 0)
                {
                    //tempArray.push_back(CreateDegenerateEdge(boundaryPoints[maxAngleIndex]));
                    tempArray.push_back(CreateDegenerateEdge(tempArray[2]->EndPoint(), tempArray[0]->StartPoint()));
                }
                else
                {
                    tempArray.insert(tempArray.begin() + maxAngleIndex, CreateDegenerateEdge(tempArray[maxAngleIndex - 1]->EndPoint(), tempArray[maxAngleIndex]->StartPoint()));
                }
                while (tempArray[3]->StartPoint().Distance(tempArray[3]->EndPoint()) > 10.0)
                {
                    tempArray.insert(tempArray.begin(), tempArray[3]);
                    tempArray.pop_back();  // 使用 pop_back 替代 erase
                }
            }
        }
        Handle(Geom_BSplineCurve) bslpineCurve1, bslpineCurve2, bslpineCurve3, bslpineCurve4;
        SurfaceModelingTool::Arrange_Coons_G0(tempArray, bslpineCurve1, bslpineCurve2, bslpineCurve3, bslpineCurve4, 10, Standard_True);
        std::vector<Handle(Geom_BSplineCurve)> aBoundarycurveArray = { bslpineCurve1 , bslpineCurve2, bslpineCurve3, bslpineCurve4 };

        // = Standard_True;
        Visualize(tempArray[0], Quantity_NOC_RED);
        Visualize(tempArray[2], Quantity_NOC_RED);
        Visualize(tempArray[1], Quantity_NOC_BLUE);
        Visualize(tempArray[3]->StartPoint(), Quantity_NOC_BLUE);
        Standard_Integer isoCount = 20;
        std::vector<Handle(Geom_BSplineCurve)> uISOcurvesArray_Initial, vISOcurvesArray_Initial;
        std::vector<Handle(Geom_BSplineCurve)> uISOcurvesArray_Final, vISOcurvesArray_Final;
        std::vector<gp_Vec> normalsOfUISOLines, normalsOfVISOLines;
        std::vector<std::vector<Standard_Real>> uKnots;
        std::vector<std::vector<Standard_Real>> vKnots;
        // 新添加代码
        std::vector<Handle(Geom_BSplineCurve)> uInternalCurve, vInternalCurve;
        Standard_Real uAngleSum = 0, vAngleSum = 0;
        Handle(Geom_BSplineSurface) aResSurface;

        //MathTool::TrimInternalCurves(anInternalBSplineCurves, aBoundarycurveArray, 20);
        // = Standard_True;
        myOccView->getContext()->RemoveAll(Standard_True);
        std::vector<Handle(Geom_BSplineCurve)> aTempCurveArray;
        aTempCurveArray.insert(aTempCurveArray.end(), aBoundarycurveArray.begin(), aBoundarycurveArray.end());
        aTempCurveArray.insert(aTempCurveArray.end(), anInternalBSplineCurves.begin(), anInternalBSplineCurves.end());
        if (MathTool::AreBSplinesCoPlanar(aTempCurveArray, 1e-3))
        {
            // 所有线共平面，直接生成Coons
            SurfaceModelingTool::Coons_G0(bslpineCurve1, bslpineCurve2, bslpineCurve3, bslpineCurve4, aResSurface);
            Visualize(aResSurface);
        }

        if (SurfaceModelingTool::GetInternalCurves(aBoundarycurveArray, anInternalBSplineCurves, uInternalCurve, vInternalCurve, uAngleSum, vAngleSum, 5))
        {
            myOccView->getContext()->RemoveAll(Standard_True);

            for (auto uCurve : uInternalCurve)
            {
                PlanarCurve u(uCurve);
                if (u.GetCurveType() == CurveType::PLANAR)
                {
                    Visualize(u.GetPlane(), Quantity_NOC_GOLD);
                }
                Visualize(u.GetCurve(), Quantity_NOC_BLUE);
            }
            for (auto vCurve : vInternalCurve)
            {
                PlanarCurve v(vCurve);
                if (v.GetCurveType() == CurveType::PLANAR)
                {
                    Visualize(v.GetPlane(), Quantity_NOC_GOLD);
                }
                Visualize(v.GetCurve(), Quantity_NOC_RED);
            }
            //Visualize(uInternalCurve, Quantity_NOC_RED);
            //Visualize(vInternalCurve, Quantity_NOC_BLUE);
            aResSurface = SurfaceModelingTool::GenerateReferSurface(aBoundarycurveArray, uInternalCurve, vInternalCurve, uAngleSum, vAngleSum, isoCount, GORDEN_TWO_DIRECTION_GORDEN);
            //Visualize(aResSurface);
        }
        myOccView->getContext()->RemoveAll(Standard_True);
        Visualize(bslpineCurve1, Quantity_NOC_RED);
        Visualize(bslpineCurve2, Quantity_NOC_RED);
        Visualize(bslpineCurve3, Quantity_NOC_RED);
        Visualize(bslpineCurve4, Quantity_NOC_RED);
        if (aResSurface.IsNull())
        {
            // 生成不了Gorden，不适用新算法，继续生成Coons
            SurfaceModelingTool::Coons_G0(bslpineCurve1, bslpineCurve2, bslpineCurve3, bslpineCurve4, aResSurface);
        }
        Visualize(aResSurface, Quantity_NOC_GOLD);
        std::vector<TopoDS_Edge> NormalArray;
        SurfaceModelingTool::GetISOCurveWithNormal(aResSurface, uISOcurvesArray_Initial, vISOcurvesArray_Initial, normalsOfUISOLines, normalsOfVISOLines, NormalArray, isoCount);

        if (MathTool::ComputeCurveCurveDistance(vISOcurvesArray_Initial[0], aBoundarycurveArray[0]) > 10.0)
        {
            std::swap(uISOcurvesArray_Initial, vISOcurvesArray_Initial);
            std::swap(normalsOfUISOLines, normalsOfVISOLines);
        }
        else if (MathTool::ComputeCurveCurveDistance(uISOcurvesArray_Initial[0], aBoundarycurveArray[1]) > 10.0)
        {
            std::swap(uISOcurvesArray_Initial, vISOcurvesArray_Initial);
            std::swap(normalsOfUISOLines, normalsOfVISOLines);
        }

        std::vector<std::pair<Handle(Geom_BSplineCurve), gp_Vec>> uInitialCurvesWithNormal;
        for (Standard_Integer m = 0; m < uISOcurvesArray_Initial.size(); m++)
        {
            uInitialCurvesWithNormal.emplace_back(uISOcurvesArray_Initial[m], normalsOfUISOLines[m]);
        }
        std::vector<std::pair<Handle(Geom_BSplineCurve), gp_Vec>> vInitialCurvesWithNormal;
        for (Standard_Integer m = 0; m < vISOcurvesArray_Initial.size(); m++)
        {
            vInitialCurvesWithNormal.emplace_back(vISOcurvesArray_Initial[m], normalsOfVISOLines[m]);
        }
        MathTool::SortBSplineCurves(uInitialCurvesWithNormal, aBoundarycurveArray[0]);
        MathTool::SortBSplineCurves(vInitialCurvesWithNormal, aBoundarycurveArray[1]);

        uISOcurvesArray_Initial.clear();
        normalsOfUISOLines.clear();
        for (std::pair<Handle(Geom_BSplineCurve), gp_Vec> uCurveWithNormal : uInitialCurvesWithNormal)
        {
            uISOcurvesArray_Initial.push_back(uCurveWithNormal.first);
            normalsOfUISOLines.emplace_back(uCurveWithNormal.second);
        }
        uISOcurvesArray_Initial.insert(uISOcurvesArray_Initial.begin(), aBoundarycurveArray[0]);
        uISOcurvesArray_Initial.push_back(aBoundarycurveArray[2]);
        normalsOfUISOLines.insert(normalsOfUISOLines.begin(), gp_Vec());
        normalsOfUISOLines.push_back(gp_Vec());

        vISOcurvesArray_Initial.clear();
        normalsOfVISOLines.clear();
        for (std::pair<Handle(Geom_BSplineCurve), gp_Vec> vCurveWithNormal : vInitialCurvesWithNormal)
        {
            vISOcurvesArray_Initial.push_back(vCurveWithNormal.first);
            normalsOfVISOLines.push_back(vCurveWithNormal.second);
        }
        vISOcurvesArray_Initial.insert(vISOcurvesArray_Initial.begin(), aBoundarycurveArray[1]);
        vISOcurvesArray_Initial.push_back(aBoundarycurveArray[3]);
        normalsOfVISOLines.insert(normalsOfVISOLines.begin(), gp_Vec());
        normalsOfVISOLines.push_back(gp_Vec());

        MathTool::ReverseIfNeeded(uISOcurvesArray_Initial);
        MathTool::ReverseIfNeeded(vISOcurvesArray_Initial);
        //构造Lofting曲面
        std::vector<TopoDS_Shape>  uLoftingSur, vLoftingSur;
        std::vector<Handle(Geom_BSplineCurve)> uLoftingCurves;
        std::vector<Handle(Geom_BSplineCurve)> vLoftingCurves;
        SurfaceModelingTool::CreateLoftingSurface(uISOcurvesArray_Initial, normalsOfUISOLines, uLoftingSur, uLoftingCurves);
        SurfaceModelingTool::CreateLoftingSurface(vISOcurvesArray_Initial, normalsOfVISOLines, vLoftingSur, vLoftingCurves);
        uLoftingSur.insert(uLoftingSur.begin(), TopoDS_Shape());
        vLoftingSur.insert(vLoftingSur.begin(), TopoDS_Shape());
        uLoftingSur.push_back(TopoDS_Shape());
        vLoftingSur.push_back(TopoDS_Shape());
        myOccView->getContext()->RemoveAll(Standard_True);
        //Visualize(NormalArray, Quantity_NOC_GOLD);
        Visualize(uISOcurvesArray_Initial, Quantity_NOC_WHITE);
        Visualize(vISOcurvesArray_Initial, Quantity_NOC_WHITE);

        Visualize(uLoftingSur);
        Visualize(vLoftingSur);
        Visualize(uLoftingCurves, Quantity_NOC_RED);
        Visualize(vLoftingCurves, Quantity_NOC_RED);
        // 生成修正的等参线
        std::vector<Handle(Geom_BSplineCurve)> uISOcurvesArray_New, vISOcurvesArray_New;
        std::vector<gp_Pnt> interPoints;
        std::vector<TopoDS_Edge> uInterpoalteTangentArray;
        std::vector<TopoDS_Edge> uInterpoalteTangentArray2;
        std::vector<TopoDS_Edge> vInterpoalteTangentArray;
        std::vector<TopoDS_Edge> vInterpoalteTangentArray2;

        std::vector<std::vector<gp_Pnt>> uInterpolatePoints;
        std::vector<std::vector<gp_Pnt>> vInterpolatePoints;
        aResSurface = 0;
        SurfaceModelingTool::LoftSurfaceIntersectWithCurve(uLoftingSur, uISOcurvesArray_Initial, anInternalBSplineCurves,
            uISOcurvesArray_New, isoCount,
            uInterpolatePoints,
            uInterpoalteTangentArray, uInterpoalteTangentArray2, aResSurface);
        myOccView->getContext()->RemoveAll(Standard_True);
        Visualize(aBoundarycurveArray, Quantity_NOC_RED);
        Visualize(uInterpoalteTangentArray, Quantity_NOC_GOLD);
        SurfaceModelingTool::LoftSurfaceIntersectWithCurve(vLoftingSur, vISOcurvesArray_Initial, anInternalBSplineCurves,
            vISOcurvesArray_New, isoCount,
            vInterpolatePoints,
            vInterpoalteTangentArray, vInterpoalteTangentArray2, aResSurface);
        myOccView->getContext()->RemoveAll(Standard_True);
        Visualize(uISOcurvesArray_New, Quantity_NOC_BLUE);
        Visualize(vISOcurvesArray_New, Quantity_NOC_RED);
        if (uInternalCurve.size() >= 4 || vInternalCurve.size() >= 4)
        {
            if (MathTool::ComputeCurveCurveDistance(vISOcurvesArray_New[0], uInternalCurve[0]) > 10 &&
                MathTool::ComputeCurveCurveDistance(vISOcurvesArray_New[0], uInternalCurve[uInternalCurve.size() - 1]) > 10)
            {
                std::swap(uISOcurvesArray_New, vISOcurvesArray_New);
            }
            if (MathTool::ComputeCurveCurveDistance(uISOcurvesArray_New[0], vInternalCurve[0]) > 10 &&
                MathTool::ComputeCurveCurveDistance(uISOcurvesArray_New[0], vInternalCurve[vInternalCurve.size() - 1]) > 10)
            {
                std::swap(uISOcurvesArray_New, vISOcurvesArray_New);
            }


            // 保留线多的方向
            if (uInternalCurve.size() > vInternalCurve.size())
            {
                uISOcurvesArray_New.clear();
                uISOcurvesArray_New = uInternalCurve;
                vISOcurvesArray_New.insert(vISOcurvesArray_New.begin(), aBoundarycurveArray[1]);
                vISOcurvesArray_New.emplace_back(aBoundarycurveArray[3]);
            }
            else
            {
                vISOcurvesArray_New.clear();
                vISOcurvesArray_New = vInternalCurve;
                uISOcurvesArray_New.insert(uISOcurvesArray_New.begin(), aBoundarycurveArray[0]);
                uISOcurvesArray_New.emplace_back(aBoundarycurveArray[2]);
            }

            MathTool::SortBSplineCurves(uISOcurvesArray_New, uISOcurvesArray_New[0]);
            MathTool::SortBSplineCurves(vISOcurvesArray_New, vISOcurvesArray_New[0]);
            MathTool::ReverseIfNeeded(uISOcurvesArray_New);
            MathTool::ReverseIfNeeded(vISOcurvesArray_New);
            myOccView->getContext()->RemoveAll(Standard_True);
            Visualize(uISOcurvesArray_New, Quantity_NOC_BLUE);
            Visualize(vISOcurvesArray_New, Quantity_NOC_RED);

            // 调用陈鑫的 Compatible
            std::for_each(uISOcurvesArray_New.begin(), uISOcurvesArray_New.end(), UniformCurve);
            std::for_each(vISOcurvesArray_New.begin(), vISOcurvesArray_New.end(), UniformCurve);
            for (auto uCurve : uInternalCurve)
            {
                PlanarCurve uPlanarCurve(uCurve);
                if (uPlanarCurve.GetCurveType() == CurveType::PLANAR)
                {
                    Visualize(uPlanarCurve.GetPlane(), Quantity_NOC_RED);
                }
                Visualize(uPlanarCurve.GetCurve(), Quantity_NOC_BLUE);
            }
            for (auto vCurve : vInternalCurve)
            {
                PlanarCurve vPlanarCurve(vCurve);
                if (vPlanarCurve.GetCurveType() == CurveType::PLANAR)
                {
                    Visualize(vPlanarCurve.GetPlane(), Quantity_NOC_RED);
                }
                Visualize(vPlanarCurve.GetCurve(), Quantity_NOC_RED);
            }
            Visualize(vISOcurvesArray_New);
            Visualize(uISOcurvesArray_New);
            CurveOperate::CompatibleWithInterPointsThree(vISOcurvesArray_New, uISOcurvesArray_New);
            Visualize(vISOcurvesArray_New, Quantity_NOC_WHITE);
            Visualize(uISOcurvesArray_New, Quantity_NOC_WHITE);
            CurveOperate::CompatibleWithInterPointsThree(uISOcurvesArray_New, vISOcurvesArray_New);
            myOccView->getContext()->RemoveAll(Standard_True);
            Visualize(vISOcurvesArray_New, Quantity_NOC_WHITE);
            Visualize(uISOcurvesArray_New, Quantity_NOC_WHITE);
            //visualize intersection point and params
            std::vector<std::vector<gp_Pnt>> pointsOnTheCurve(vISOcurvesArray_New.size());
            std::vector<std::vector<Standard_Real>> paramsOnTheCurve(vISOcurvesArray_New.size());
            for (Standard_Integer i = 0; i < vISOcurvesArray_New.size(); i++) {
                if (CurveOperate::isDegenerate(vISOcurvesArray_New[i])) {
                    pointsOnTheCurve[i] = { vISOcurvesArray_New[i]->Pole(1) };
                    paramsOnTheCurve[i] = { vISOcurvesArray_New[i]->FirstParameter() };
                    continue;
                }
                std::tie(pointsOnTheCurve[i], paramsOnTheCurve[i]) = CurveOperate::CalCurvesInterPointsParamsToCurve(uISOcurvesArray_New, vISOcurvesArray_New[i]);
            }
            // 创建一个存储pair的vector
            std::vector<std::pair<gp_Pnt, Standard_Real>> pointParamPairs;
            for (Standard_Integer i = 0; i < pointsOnTheCurve.size(); i++) {
                for (Standard_Integer j = 0; j < pointsOnTheCurve[i].size(); j++) {
                    pointParamPairs.emplace_back(std::make_pair(pointsOnTheCurve[i][j], paramsOnTheCurve[i][j]));
                }
            }
            Visualize(pointParamPairs, Quantity_NOC_GOLD);
            //visualize intersection point and params
            std::vector<std::vector<gp_Pnt>> pointsOnTheCurve1(uISOcurvesArray_New.size());
            std::vector<std::vector<Standard_Real>> paramsOnTheCurve1(uISOcurvesArray_New.size());
            for (Standard_Integer i = 0; i < uISOcurvesArray_New.size(); i++) {
                if (CurveOperate::isDegenerate(uISOcurvesArray_New[i])) {
                    pointsOnTheCurve1[i] = { uISOcurvesArray_New[i]->Pole(1) };
                    paramsOnTheCurve1[i] = { uISOcurvesArray_New[i]->FirstParameter() };
                    continue;
                }
                std::tie(pointsOnTheCurve1[i], paramsOnTheCurve1[i]) = CurveOperate::CalCurvesInterPointsParamsToCurve(vISOcurvesArray_New, uISOcurvesArray_New[i]);
            }
            // 创建一个存储pair的vector
            std::vector<std::pair<gp_Pnt, Standard_Real>> pointParamPairs1;
            for (Standard_Integer i = 0; i < pointsOnTheCurve1.size(); i++) {
                for (Standard_Integer j = 0; j < pointsOnTheCurve1[i].size(); j++) {
                    pointParamPairs1.emplace_back(std::make_pair(pointsOnTheCurve1[i][j], paramsOnTheCurve1[i][j]));
                }
            }
            myOccView->getContext()->RemoveAll(Standard_True);
            Visualize(uISOcurvesArray_New, Quantity_NOC_BLUE);
            Visualize(vISOcurvesArray_New, Quantity_NOC_RED);
            Visualize(pointParamPairs1, Quantity_NOC_GOLD);
            // 调用胡新宇的 Gorden
            TopoDS_Face aGordenFace;
            std::vector<gp_Pnt> upoints, vpoints;
            std::vector<Standard_Real> uparams, vparams;
            Standard_Integer vIndex, uIndex;
            vIndex = CurveOperate::isDegenerate(vISOcurvesArray_New.front()) ? vISOcurvesArray_New.size() - 1 : 0;
            uIndex = CurveOperate::isDegenerate(uISOcurvesArray_New.front()) ? uISOcurvesArray_New.size() - 1 : 0;
            std::tie(upoints, uparams) = CurveOperate::CalCurvesInterPointsParamsToCurve(uISOcurvesArray_New, vISOcurvesArray_New[vIndex]);
            std::tie(vpoints, vparams) = CurveOperate::CalCurvesInterPointsParamsToCurve(vISOcurvesArray_New, uISOcurvesArray_New[uIndex]);
            // 排序操作
            std::vector<std::pair<Standard_Real, Handle(Geom_BSplineCurve)>> combinedv;
            for (size_t i = 0; i < vparams.size(); ++i) {
                combinedv.emplace_back(vparams[i], vISOcurvesArray_New[i]);
            }
            std::sort(combinedv.begin(), combinedv.end(), [](const auto& a, const auto& b) {
                return a.first < b.first;
                });
            for (size_t i = 0; i < combinedv.size(); ++i) {
                vparams[i] = combinedv[i].first;
                vISOcurvesArray_New[i] = combinedv[i].second;
            }

            std::vector<std::pair<Standard_Real, Handle(Geom_BSplineCurve)>> combinedu;
            for (size_t i = 0; i < uparams.size(); ++i) {
                combinedu.emplace_back(uparams[i], uISOcurvesArray_New[i]);
            }
            std::sort(combinedu.begin(), combinedu.end(), [](const auto& a, const auto& b) {
                return a.first < b.first;
                });
            for (size_t i = 0; i < combinedu.size(); ++i) {
                uparams[i] = combinedu[i].first;
                uISOcurvesArray_New[i] = combinedu[i].second;
            }
            GordenSurface::BuildMyGordonSurf(uISOcurvesArray_New, vISOcurvesArray_New, uparams, vparams, aGordenFace);
            Handle(Geom_Surface) aGeomSurface = BRep_Tool::Surface(aGordenFace);
            Handle(Geom_BSplineSurface) aBSplineSurface = Handle(Geom_BSplineSurface)::DownCast(aGeomSurface);
            Visualize(aBSplineSurface);
            // 直接不进行后面的操作，计算下一个case
            continue;;
        }
        myOccView->getContext()->RemoveAll(Standard_True);
        // 根据u、v等参线之间的交点，生成最终等参线
        interPoints.clear();
        std::vector<gp_Pnt> boundaryPoints;
        std::vector<Handle(Geom_BSplineSurface)> surfaceArray;
        std::vector<TopoDS_Edge> TangentArray;
        std::string gordenSurf1 = filename + "gordonSurf1.step";
        std::string gordenSurf2 = filename + "gordonSurf2.step";
        std::string gordenSurf3 = filename + "gordonSurf3.step";
        std::string gordenSurf4 = filename + "gordonSurf4.step";
        SurfaceModelingTool::LoadBSplineSurfaces(gordenSurf1, surfaceArray);
        SurfaceModelingTool::LoadBSplineSurfaces(gordenSurf2, surfaceArray);
        SurfaceModelingTool::LoadBSplineSurfaces(gordenSurf3, surfaceArray);
        SurfaceModelingTool::LoadBSplineSurfaces(gordenSurf4, surfaceArray);
        Visualize(surfaceArray);
        SurfaceModelingTool::CreateFinalISOCurves(uISOcurvesArray_New, vISOcurvesArray_New,
            uISOcurvesArray_Final, vISOcurvesArray_Final,
            uInterpolatePoints, vInterpolatePoints,
            uKnots, vKnots,
            boundaryPoints, interPoints,
            isoCount,
            TangentArray, surfaceArray);
        Visualize(TangentArray, Quantity_NOC_GOLD);
        uISOcurvesArray_Final.insert(uISOcurvesArray_Final.begin(), aBoundarycurveArray[0]);
        uISOcurvesArray_Final.push_back(aBoundarycurveArray[2]);
        vISOcurvesArray_Final.insert(vISOcurvesArray_Final.begin(), aBoundarycurveArray[1]);
        vISOcurvesArray_Final.push_back(aBoundarycurveArray[3]);
        MathTool::SortBSplineCurves(uISOcurvesArray_Final, uISOcurvesArray_Final[0]);
        MathTool::SortBSplineCurves(vISOcurvesArray_Final, vISOcurvesArray_Final[0]);
        MathTool::ReverseIfNeeded(uISOcurvesArray_Final);
        MathTool::ReverseIfNeeded(vISOcurvesArray_Final);
        for (Standard_Integer i = 0; i < uISOcurvesArray_Final.size(); i++)
        {
            Visualize(std::make_pair(uISOcurvesArray_Final[i]->Pole(1), uISOcurvesArray_Final[i]->FirstParameter()));
        }
        for (auto uCurves : uISOcurvesArray_Final)
        {
            Visualize(uCurves->StartPoint(), Quantity_NOC_GOLD);
            Visualize(uCurves->EndPoint(), Quantity_NOC_AZURE);
        }

        for (auto vCurves : vISOcurvesArray_Final)
        {
            Visualize(vCurves->StartPoint(), Quantity_NOC_GOLD);
            Visualize(vCurves->EndPoint(), Quantity_NOC_AZURE);
        }
        Visualize(uISOcurvesArray_Final, Quantity_NOC_RED);
        Visualize(vISOcurvesArray_Final, Quantity_NOC_BLUE);
        // 调用胡新宇的 Gorden
        CurveOperate::CompatibleWithInterPointsThree(vISOcurvesArray_Final, uISOcurvesArray_Final);
        CurveOperate::CompatibleWithInterPointsThree(uISOcurvesArray_Final, vISOcurvesArray_Final);
        TopoDS_Face aGordenFace;
        std::vector<gp_Pnt> upoints, vpoints;
        std::vector<Standard_Real> uparams, vparams;
        Standard_Integer vIndex, uIndex;
        vIndex = CurveOperate::isDegenerate(vISOcurvesArray_Final.front(), 10) ? vISOcurvesArray_Final.size() - 1 : 0;
        uIndex = CurveOperate::isDegenerate(uISOcurvesArray_Final.front(), 10) ? uISOcurvesArray_Final.size() - 1 : 0;
        LogPrint::PrintInfo(vISOcurvesArray_Final);
        std::cout << "---------------------------------------------------" << std::endl;
        LogPrint::PrintInfo(uISOcurvesArray_Final);
        std::tie(upoints, uparams) = CurveOperate::CalCurvesInterPointsParamsToCurve(uISOcurvesArray_Final, vISOcurvesArray_Final[vIndex]);
        std::tie(vpoints, vparams) = CurveOperate::CalCurvesInterPointsParamsToCurve(vISOcurvesArray_Final, uISOcurvesArray_New[uIndex]);
        // 排序操作
        std::vector<std::pair<Standard_Real, Handle(Geom_BSplineCurve)>> combinedv;
        for (size_t i = 0; i < vparams.size(); ++i) {
            combinedv.emplace_back(vparams[i], vISOcurvesArray_Final[i]);
        }
        std::sort(combinedv.begin(), combinedv.end(), [](const auto& a, const auto& b) {
            return a.first < b.first;
            });
        for (size_t i = 0; i < combinedv.size(); ++i) {
            vparams[i] = combinedv[i].first;
            vISOcurvesArray_Final[i] = combinedv[i].second;
        }

        std::vector<std::pair<Standard_Real, Handle(Geom_BSplineCurve)>> combinedu;
        for (size_t i = 0; i < uparams.size(); ++i) {
            combinedu.emplace_back(uparams[i], uISOcurvesArray_Final[i]);
        }
        std::sort(combinedu.begin(), combinedu.end(), [](const auto& a, const auto& b) {
            return a.first < b.first;
            });
        for (size_t i = 0; i < combinedu.size(); ++i) {
            uparams[i] = combinedu[i].first;
            uISOcurvesArray_Final[i] = combinedu[i].second;
        }
        GordenSurface::BuildMyGordonSurf(uISOcurvesArray_Final, vISOcurvesArray_Final, uparams, vparams, aGordenFace);
        Handle(Geom_Surface) aGeomSurface = BRep_Tool::Surface(aGordenFace);
        Handle(Geom_BSplineSurface) aBSplineSurface = Handle(Geom_BSplineSurface)::DownCast(aGeomSurface);
        Visualize(aBSplineSurface);
        std::string uIsoExportFilename = "D:/GordenModels/NewModels/Curves/";
        uIsoExportFilename += "coons_";
        uIsoExportFilename += std::to_string(i);
        uIsoExportFilename += "_uIsoCurves.brep";
        std::string vIsoExportFilename = "D:/GordenModels/NewModels/Curves/";
        vIsoExportFilename += "coons_";
        vIsoExportFilename += std::to_string(i);
        vIsoExportFilename += "_vIsoCurves.brep";

        SurfaceModelingTool::ExportBSplineCurves(uISOcurvesArray_Final, uIsoExportFilename);
        SurfaceModelingTool::ExportBSplineCurves(vISOcurvesArray_Final, vIsoExportFilename);
        for (auto boundaryPoint : boundaryPoints)
        {
            interPoints.push_back(boundaryPoint);
        }
        interPoints.push_back(uISOcurvesArray_Final[0]->StartPoint());
        interPoints.push_back(uISOcurvesArray_Final[0]->EndPoint());
        interPoints.push_back(uISOcurvesArray_Final[uISOcurvesArray_Final.size() - 1]->StartPoint());
        interPoints.push_back(uISOcurvesArray_Final[uISOcurvesArray_Final.size() - 1]->EndPoint());

        //Visualize(interPoints, Quantity_NOC_GOLD);
        // 遍历 u(v)ISOcurvesArray_Final 进行可视化
        Visualize(uISOcurvesArray_Final);
        Visualize(vISOcurvesArray_Final);
        auto ExportPointsToBREP = [](const std::vector<gp_Pnt>& boundaryPoints, const std::string& filename)
            {
                // 创建一个复合体以存储所有顶点
                TopoDS_Compound compound;
                BRep_Builder builder;
                builder.MakeCompound(compound);

                // 将 gp_Pnt 转换为 TopoDS_Vertex 并添加到复合体
                for (const auto& point : boundaryPoints)
                {
                    TopoDS_Vertex vertex = BRepBuilderAPI_MakeVertex(point);
                    builder.Add(compound, vertex);
                }

                // 将复合体保存到 BREP 文件
                if (BRepTools::Write(compound, filename.c_str()))
                {
                    std::cout << "成功导出到 BREP 文件: " << filename << std::endl;
                }
                else {
                    std::cerr << "导出 BREP 文件失败！" << std::endl;
                }
            };

        ExportPointsToBREP(interPoints, filename + std::string("points.brep"));

        TopoDS_Compound UResult;
        BRep_Builder builder1;
        builder1.MakeCompound(UResult);
        for (auto curve : uISOcurvesArray_Final)
        {
            TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(curve);
            builder1.Add(UResult, edge);
        }
        TopoDS_Compound VResult;
        BRep_Builder builder;
        builder.MakeCompound(VResult);
        for (auto curve : vISOcurvesArray_Final)
        {
            TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(curve);
            builder.Add(VResult, edge);
        }

        std::string UresultPath = filename + "UResult.brep"; std::string VresultPath = filename + "VResult.brep";
        BRepTools::Write(UResult, UresultPath.c_str()); BRepTools::Write(VResult, VresultPath.c_str());

        SurfaceModelingTool tool;
        std::string knotsPath = filename + "knots.txt";
        tool.setKnotsOutputPath(knotsPath.c_str());
        // 检查文件是否存在，如果存在，清空文件内容
        std::ifstream checkFile(tool.getKnotsOuputPath());
        if (checkFile.is_open())
        {
            // 关闭检查文件的输入流
            checkFile.close();
            // 清空文件内容，覆盖写
            std::ofstream clearFile(tool.getKnotsOuputPath(), std::ios::trunc);
            clearFile.close();
        }

        tool.ContextToTxt("U:");
        for (auto debugKnots : uKnots)
            tool.KnotsToTxt(debugKnots);

        tool.ContextToTxt("------------------------------\nV:");
        for (auto debugKnots : vKnots)
            tool.KnotsToTxt(debugKnots);
    }
}
