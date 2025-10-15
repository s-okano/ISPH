#include <vtkSmartPointer.h>
#include <vtkIntArray.h>
#include <vtkFloatArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPointData.h>
#include <vtkVertex.h>
#include <vtkCellArray.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>

vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid;
extern "C" {

  std::string& ltrim(std::string& str, const std::string& chars = "\t\n\v\f\r ")
  {
      str.erase(0, str.find_first_not_of(chars));
      return str;
  }
   
  std::string& rtrim(std::string& str, const std::string& chars = "\t\n\v\f\r ")
  {
      str.erase(str.find_last_not_of(chars) + 1);
      return str;
  }
   
  std::string& alltrim(std::string& str, const std::string& chars = "\t\n\v\f\r ")
  {
      return ltrim(rtrim(str, chars), chars);
  }

  void initialize()
  {
    unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  }

  void setPoints(int np, float *coord)
  {
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for(int i = 0 ; i < 3*np ; i+=3) points->InsertNextPoint(coord[i], coord[i+1], coord[i+2]);
    unstructuredGrid->SetPoints(points);
    
    vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New();
    
    for(int i = 0 ; i < np ; i++)
    {
        vtkSmartPointer<vtkVertex> element = vtkSmartPointer<vtkVertex>::New();
        element->GetPointIds()->SetId(0,  i);
        cellArray->InsertNextCell(element);
    }
    int *cellTypes = new int[np];
    for(int i = 0 ; i < np ; i++) cellTypes[i] = VTK_VERTEX;

    unstructuredGrid->SetCells(cellTypes, cellArray);
  }


  void setPointFloatScalarValue(int np, float *value, int len_name, const char* name)
  {
    vtkSmartPointer<vtkFloatArray> pointData = vtkSmartPointer<vtkFloatArray>::New();
    std::string str(name,0,len_name); // ISO_C_BINDING couldn't send correct char-array size.
    pointData->SetName(str.c_str());
    pointData->SetNumberOfComponents(1);
    pointData->SetNumberOfValues(np);
    for(int i = 0; i < np; i++) pointData->SetValue(i, value[i]);
    unstructuredGrid->GetPointData()->AddArray(pointData);
  }

  void setPointIntScalarValue(int np, int *value, int len_name, const char* name)
  {
    vtkSmartPointer<vtkIntArray> pointData = vtkSmartPointer<vtkIntArray>::New();
    std::string str(name,0,len_name);
    pointData->SetName(str.c_str());
    pointData->SetNumberOfComponents(1);
    pointData->SetNumberOfValues(np);
    for(int i = 0; i < pointData->GetNumberOfTuples(); i++) pointData->SetValue(i, value[i]);
    unstructuredGrid->GetPointData()->AddArray(pointData);
  }

  
  void setPointFloatVectorValue(int np, float *value, int len_name, const char* name)
  {
    vtkSmartPointer<vtkFloatArray> pointData = vtkSmartPointer<vtkFloatArray>::New();
    std::string str(name,0,len_name);
    pointData->SetName(str.c_str());
    pointData->SetNumberOfComponents(3);
    pointData->SetNumberOfValues(3*np);
    for(int i = 0; i < np; i++) pointData->SetTuple3(i, value[3*i],value[3*i+1],value[3*i+2]);
    unstructuredGrid->GetPointData()->AddArray(pointData);
  }

  void output(int len_filename, const char* filename)
  {
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    //writer->SetDataModeToAscii(); // DEBUG : ASCIIモード
    std::string str(filename,0,len_filename);
    writer->SetFileName(str.c_str());
#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(unstructuredGrid);
#else
    writer->SetInputData(unstructuredGrid);
#endif
    writer->Write();
  }


  void finalize()
  {
    unstructuredGrid = NULL;
  }
}
