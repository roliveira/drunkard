#ifndef IO_WRITER_H_
#define IO_WRITER_H_


#include <iostream>
#include <vector>

#include <vtkSmartPointer.h>
#include <vtkProperty.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkImageData.h>

#include "xtensor/xview.hpp"
#include "xtensor/xarray.hpp"
#include "xtensor/xadapt.hpp"
#include "xtensor/xbuilder.hpp"

#include "utils/index.hpp"


template<typename T>
inline void ToRawFile(std::vector<T> data, const char* fname) {

    std::ofstream f(fname, std::ios::binary | std::ios::out);

    if (f.is_open()) {
        for (T i : data) f.write(reinterpret_cast<char*>(&i), sizeof(T));
        f.close();
    }

}

void ToVTI(const char *fname, std::vector<double> data, std::vector<int> shape, std::vector<double> spacing, std::vector<double> origin) {

    std::vector<double> s {
        static_cast<double>(spacing[0]),
        static_cast<double>(spacing[1]),
        static_cast<double>(spacing[2])
    };
    
    std::vector<double> o {
        static_cast<double>(origin[0]),
        static_cast<double>(origin[1]),
        static_cast<double>(origin[2])
    };

    vtkSmartPointer<vtkImageData> image_data =
        vtkSmartPointer<vtkImageData>::New();

    image_data->SetDimensions(shape.data());
    image_data->SetSpacing(s.data());
    image_data->SetOrigin(o.data());
    image_data->AllocateScalars(VTK_DOUBLE, 1);

    for (int z = 0; z < shape[2]; ++z) {
        for (int y = 0; y < shape[1]; ++y) {
            for (int x = 0; x < shape[0]; ++x) {
                int linear_index = sub2ind({x, y, z}, shape);
                double* voxel = static_cast<double*>(image_data->GetScalarPointer(x, y, z));
                voxel[0] = data[linear_index];
            }
        }
    }
    
    vtkSmartPointer<vtkXMLImageDataWriter> writer =
        vtkSmartPointer<vtkXMLImageDataWriter>::New();

    writer->SetFileName(fname);
    writer->SetInputData(image_data);
    writer->Write();

}

void ToVTI(const char *fname, xt::xarray<double> data, std::vector<int> shape, std::vector<double> spacing, std::vector<double> origin)
{
    std::vector<double> s {
        static_cast<double>(spacing[0]),
        static_cast<double>(spacing[1]),
        static_cast<double>(spacing[2])
    };
    
    std::vector<double> o {
        static_cast<double>(origin[0]),
        static_cast<double>(origin[1]),
        static_cast<double>(origin[2])
    };

    vtkSmartPointer<vtkImageData> image_data =
        vtkSmartPointer<vtkImageData>::New();

    image_data->SetDimensions(shape.data());
    image_data->SetSpacing(s.data());
    image_data->SetOrigin(o.data());
    image_data->AllocateScalars(VTK_DOUBLE, 1);

    for (int iz = 0; iz < shape[2]; ++iz) {
        for (int iy = 0; iy < shape[1]; ++iy) {
            for (int ix = 0; ix < shape[0]; ++ix) {
                image_data->SetScalarComponentFromDouble(ix, iy, iz, 0, data(ix, iy, iz));
                // double* voxel = static_cast<double*>(image_data->GetScalarPointer(ix, iy, iz));
                // voxel[0] = data(ix, iy, iz);
            }
        }
    }
    
    vtkSmartPointer<vtkXMLImageDataWriter> writer =
        vtkSmartPointer<vtkXMLImageDataWriter>::New();

    writer->SetFileName(fname);
    writer->SetInputData(image_data);
    writer->Write();
}

void Vector2VTI(const char *fname, xt::xarray<double> data, std::vector<int> shape, std::vector<double> spacing, std::vector<double> origin) {

    int dir = shape.size();

    std::vector<double> s {
        static_cast<double>(spacing[0]),
        static_cast<double>(spacing[1]),
        static_cast<double>(spacing[2])
    };
    
    std::vector<double> o {
        static_cast<double>(origin[0]),
        static_cast<double>(origin[1]),
        static_cast<double>(origin[2])
    };

    vtkSmartPointer<vtkImageData> image_data =
        vtkSmartPointer<vtkImageData>::New();
    image_data->SetDimensions(shape.data());
    image_data->SetSpacing(s.data());
    image_data->SetOrigin(o.data());
    image_data->AllocateScalars(VTK_DOUBLE, dir);
    
    for (int iz = 0; iz < shape[2]; ++iz) 
    {
        for (int iy = 0; iy < shape[1]; ++iy) 
        {
            for (int ix = 0; ix < shape[0]; ++ix) 
            {
                double* voxel = static_cast<double*>(image_data->GetScalarPointer(ix, iy, iz));

                for (int d=0; d<dir; ++d) 
                {
                    voxel[d] = data(ix, iy, iz, d);
                }
            }
        }
    }
    
    vtkSmartPointer<vtkXMLImageDataWriter> writer =
        vtkSmartPointer<vtkXMLImageDataWriter>::New();
    writer->SetFileName(fname);
    writer->SetInputData(image_data);
    writer->Write();
}

void Vector2VTI(const char *fname, std::vector<double> data, std::vector<int> shape, std::vector<double> spacing, std::vector<double> origin) {

    int dim = shape.size();

    std::vector<double> s {
        static_cast<double>(spacing[0]),
        static_cast<double>(spacing[1]),
        static_cast<double>(spacing[2])
    };
    
    std::vector<double> o {
        static_cast<double>(origin[0]),
        static_cast<double>(origin[1]),
        static_cast<double>(origin[2])
    };

    vtkSmartPointer<vtkImageData> image_data =
        vtkSmartPointer<vtkImageData>::New();
    image_data->SetDimensions(shape.data());
    image_data->SetSpacing(s.data());
    image_data->SetOrigin(o.data());
    image_data->AllocateScalars(VTK_DOUBLE, dim);
    
    for (int z = 0; z < shape[2]; ++z) {
        for (int y = 0; y < shape[1]; ++y) {
            for (int x = 0; x < shape[0]; ++x) {
                int linear_index = sub2ind({x, y, z}, shape);
                double* voxel = static_cast<double*>(image_data->GetScalarPointer(x, y, z));
                for (int d=0; d<dim; ++d) voxel[d] = data[linear_index*dim + d];
            }
        }
    }
    
    vtkSmartPointer<vtkXMLImageDataWriter> writer =
        vtkSmartPointer<vtkXMLImageDataWriter>::New();
    writer->SetFileName(fname);
    writer->SetInputData(image_data);
    writer->Write();
}

void Tensor2VTI(const char *fname, std::vector<double> data, std::vector<int> shape, std::vector<double> spacing, std::vector<double> origin) {

    int dim = 2 * shape.size();

    std::vector<double> s {
        static_cast<double>(spacing[0]),
        static_cast<double>(spacing[1]),
        static_cast<double>(spacing[2])
    };
    
    std::vector<double> o {
        static_cast<double>(origin[0]),
        static_cast<double>(origin[1]),
        static_cast<double>(origin[2])
    };

    vtkSmartPointer<vtkImageData> image_data =
        vtkSmartPointer<vtkImageData>::New();

    image_data->SetDimensions(shape.data());
    image_data->SetSpacing(s.data());
    image_data->SetOrigin(o.data());
    image_data->AllocateScalars(VTK_DOUBLE, dim);
    
    for (int z = 0; z < shape[2]; ++z) {
        for (int y = 0; y < shape[1]; ++y) {
            for (int x = 0; x < shape[0]; ++x) {
                int linear_index = sub2ind({x, y, z}, shape);
                double* voxel = static_cast<double*>(image_data->GetScalarPointer(x, y, z));
                for (int d=0; d<dim; ++d) voxel[d] = data[linear_index*dim + d];
            }
        }
    }
    
    vtkSmartPointer<vtkXMLImageDataWriter> writer =
        vtkSmartPointer<vtkXMLImageDataWriter>::New();

    writer->SetFileName(fname);
    writer->SetInputData(image_data);
    writer->Write();

}

void Tensor2VTI(const char *fname, xt::xarray<double> data, std::vector<int> shape, std::vector<double> spacing, std::vector<double> origin) {

    int dir = 3;

    std::vector<double> s {
        static_cast<double>(spacing[0]),
        static_cast<double>(spacing[1]),
        static_cast<double>(spacing[2])
    };
    
    std::vector<double> o {
        static_cast<double>(origin[0]),
        static_cast<double>(origin[1]),
        static_cast<double>(origin[2])
    };

    vtkSmartPointer<vtkImageData> image_data =
        vtkSmartPointer<vtkImageData>::New();

    image_data->SetDimensions(shape.data());
    image_data->SetSpacing(s.data());
    image_data->SetOrigin(o.data());
    image_data->AllocateScalars(VTK_DOUBLE, dir);

    for (int iz = 0; iz < shape[2]; ++iz) 
    {
        for (int iy = 0; iy < shape[1]; ++iy) 
        {
            for (int ix = 0; ix < shape[0]; ++ix) 
            {
                for (int d=0; d<dir; ++d) 
                {
                    image_data->SetScalarComponentFromDouble(ix, iy, iz, d, data(ix, iy, iz, 2*d+1));
                    // double* voxel = static_cast<double*>(image_data->GetScalarPointer(ix, iy, iz));
                    // voxel[d] = data(ix, iy, iz, 2*d);
                }
            }
        }
    }
    
    vtkSmartPointer<vtkXMLImageDataWriter> writer =
        vtkSmartPointer<vtkXMLImageDataWriter>::New();

    writer->SetFileName(fname);
    writer->SetInputData(image_data);
    writer->Write();

}


#endif  // IO_WRITER_H_
