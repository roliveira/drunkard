#ifndef IO_READER_H_
#define IO_READER_H_


#include <iostream>
#include <vector>
#include <iterator>

#include <vtkSmartPointer.h>
#include <vtkProperty.h>
#include <vtkXMLImageDataReader.h>
#include <vtkImageData.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>

#include "xtensor/xview.hpp"
#include "xtensor/xarray.hpp"
#include "xtensor/xadapt.hpp"
#include "xtensor/xbuilder.hpp"
#include "xtensor/xnpy.hpp"

#include "utils/index.hpp"


template<typename T>
inline std::vector<T> FromRawFile(const char* fname, int linear_shape) {
    std::vector<T> out;
    out.reserve(linear_shape);
    
    std::ifstream ifs(fname, std::ios::binary);

    if (ifs.is_open()) {
        std::copy(std::istream_iterator<char>(ifs), std::istream_iterator<char>(), std::back_inserter(out));
        ifs.close();
    }

    return out;
}

inline std::vector<int> FromTextFile(const char *fname) {
    std::vector<int> out;
    std::ifstream ifs(fname);

    if (ifs.is_open()) {
        std::copy(std::istream_iterator<char>(ifs), std::istream_iterator<char>(), std::back_inserter(out));
        ifs.close();
    }

    return out;
}

/*
inline std::vector<int> FromVTI(const char *fname) {
    vtkSmartPointer<vtkXMLImageDataReader> reader =
        vtkSmartPointer<vtkXMLImageDataReader>::New();
    reader->SetFileName(fname);
    reader->Update();

    vtkSmartPointer<vtkImageData> image_data =
        vtkSmartPointer<vtkImageData>::New();
    image_data->ShallowCopy(reader->GetOutput());
    
    double *range = image_data->GetScalarRange();
    int *dims = image_data->GetDimensions();
    std::vector<int> shape(dims, dims + 3);
          
    std::vector<int> data;
    data.resize(shape[0] * shape[1] * shape[2]);

    for (int z = 0; z < shape[2]; ++z) {
        for (int y = 0; y < shape[1]; ++y) {
            for (int x = 0; x < shape[0]; ++x) {
                int linear_index = sub2ind({x, y, z}, std::vector<int>(shape));
                int* voxel = static_cast<int*>(image_data->GetScalarPointer(x, y, z));
                if   (voxel[0] >= range[1]) data[linear_index] = 1;
                else                        data[linear_index] = 0;
            }
        }
    }

    return data;
}
*/

inline std::vector<float> FromVTI(const char *fname) {
    vtkSmartPointer<vtkXMLImageDataReader> reader =
        vtkSmartPointer<vtkXMLImageDataReader>::New();
    reader->SetFileName(fname);
    reader->Update();

    vtkSmartPointer<vtkImageData> image_data =
        vtkSmartPointer<vtkImageData>::New();
    image_data->DeepCopy(reader->GetOutput());

    int *dims = image_data->GetDimensions();
    std::vector<int> shape(dims, dims + 3);
          
    std::vector<float> data;
    data.resize(shape[0] * shape[1] * shape[2]);
    
    int n = image_data->GetCellData()->GetArray(0)->GetNumberOfTuples();

    for (int i=0; i<n; ++i)
    {
        data[i] = *image_data->GetCellData()->GetArray(0)->GetTuple(i);
    }

    return data;
}

xt::xarray<float> FromVTIToXArray(const char *fname) {
    vtkSmartPointer<vtkXMLImageDataReader> reader =
        vtkSmartPointer<vtkXMLImageDataReader>::New();
    reader->SetFileName(fname);
    reader->Update();

    vtkSmartPointer<vtkImageData> image_data =
        vtkSmartPointer<vtkImageData>::New();
    image_data->DeepCopy(reader->GetOutput());

    int *dims = image_data->GetDimensions();
    std::vector<int> shape(dims, dims+3);

    xt::xarray<float> data = xt::zeros<float>(shape);
    
    for (int iz=0; iz<shape[2]; ++iz)
    {
        for (int iy=0; iy<shape[1]; ++iy)
        {
            for (int ix=0; ix<shape[0]; ++ix)
            {
                data(ix, iy, iz) = image_data->GetScalarComponentAsFloat(ix, iy, iz, 0);
            }    
        }   
    }

    return data;
}

// std::vector<int> FromTextFile(const char *fname) {

//     vtkSmartPointer<vtkDelimitedTextReader> reader =
//         vtkSmartPointer<vtkDelimitedTextReader>::New();
//     reader->SetFileName(fname);
//     reader->DetectNumericColumnsOn();
//     reader->SetFieldDelimiterCharacters(" ");
//     reader->Update();

//     reader->GetOutput();

//     vtkTable* table = reader->GetOutput();

//     int nrows = table->GetNumberOfRows();
//     int ncols = table->GetNumberOfColumns();

//     int k = 0;
//     std::vector<int> out;
//     out.reserve(nrows * ncols);

//     for (vtkIdType i = 0; i < nrows; ++i) {
//         for (vtkIdType j = 0; j < ncols; ++j) {
//             out.push_back((table->GetValue(i, j)).ToInt());
//             ++k;
//         }
//     }

//     return out;
// }


#endif  // IO_READER_H_
