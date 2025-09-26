#include </usr/include/hdf5/serial/hdf5.h>
#include "hdf5.h"
#include "io.h"

#ifdef HDF5
void Write_hdf5 ( bool verbose )
{
    hid_t hdf5_file = 0, hdf5_grp[6], hdf5_headergrp = 0;
    hid_t hdf5_datatype = 0, hdf5_dataspace = 0, hdf5_dataset = 0;
    hsize_t dims[2];
    int rank = 0;
    char buf[512];
    int type;
    
    /* Set Header */
    Header.npart[0] = Param.Npart;
    Header.mass[0] = Problem.Mpart;
    Header.npartTotal[0] = Header.npart[0];
    Header.time = 0.0;
    Header.redshift = 0.0;

    for ( int i = 1; i < 6; i++ ) {
        Header.npart[i] = Header.mass[i] = Header.npartTotal[i] = 0;
    }

    Header.num_files = 1;
    Header.BoxSize = Problem.Boxsize[0];

    if ( verbose == 1 ) {
        printf ( "Output : \n"
                 "   File Name = %s.hdf5\n"
                 , Problem.Name);
    }


    //sprintf(buf, strcat(Problem.Name, ".hdf5"));
    hdf5_file = H5Fcreate(strcat(Problem.Name, ".hdf5"), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    hdf5_headergrp = H5Gcreate(hdf5_file, "/Header", 0);

    write_header_attributes_in_hdf5(hdf5_headergrp);

    for(type = 0; type < 6; type++)
    {
        if(Header.npart[type] > 0)
        {
            sprintf(buf, "/PartType%d", type);
            hdf5_grp[type] = H5Gcreate(hdf5_file, buf, 0);
        }
    }


	for(type = 0; type < 6; type++)
	{
        if(Header.npart[type] > 0)
        {  
            dims[0] = Header.npart[type];
            
            strcpy(buf, "InternalEnergy");
            dims[1] = 1;
            rank=1;
            
            float* u_buffer;
            u_buffer = (float*) malloc ( dims[0] * dims[1] * sizeof(float) );

            for (int i = 0; i < dims[0]; i++ )
                u_buffer[i] = (float) SphP[i].U;

            hdf5_dataspace = H5Screate_simple(rank, dims, NULL);
            hdf5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);
            hdf5_dataset = H5Dcreate(hdf5_grp[type], buf, hdf5_datatype,
                    hdf5_dataspace, H5P_DEFAULT);
            H5Dwrite(hdf5_dataset, hdf5_datatype, H5S_ALL, H5S_ALL,
                    H5P_DEFAULT, u_buffer);
            H5Dclose(hdf5_dataset);
            H5Tclose(hdf5_datatype);
            H5Sclose(hdf5_dataspace);

            strcpy(buf, "Density");
            dims[1] = 1;
            rank=1;
            
            float* d_buffer;
            d_buffer = (float*) malloc ( dims[0] * dims[1] * sizeof(float) );

            for (int i = 0; i < dims[0]; i++ )
                d_buffer[i] = (float) SphP[i].Rho;

            hdf5_dataspace = H5Screate_simple(rank, dims, NULL);
            hdf5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);
            hdf5_dataset = H5Dcreate(hdf5_grp[type], buf, hdf5_datatype,
                    hdf5_dataspace, H5P_DEFAULT);
            H5Dwrite(hdf5_dataset, hdf5_datatype, H5S_ALL, H5S_ALL,
                    H5P_DEFAULT, d_buffer);
            H5Dclose(hdf5_dataset);
            H5Tclose(hdf5_datatype);
            H5Sclose(hdf5_dataspace);

            strcpy(buf, "Velocities");
            rank=2;
            dims[1]=3;

            float* vel_buffer;
            vel_buffer = (float*) malloc ( dims[0] * dims[1] * sizeof(float) );

            for (int i = 0; i < dims[0]; i++ )
                for (int j = 0; j < dims[1]; j++ )
                {
                    int ipart = dims[1] * i + j;
                    vel_buffer[ipart] = (float) P[i].Vel[j];
                }
            
            hdf5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);
            hdf5_dataspace = H5Screate_simple(rank, dims, NULL);
            hdf5_dataset = H5Dcreate(hdf5_grp[type], buf, hdf5_datatype,
                    hdf5_dataspace, H5P_DEFAULT);
            H5Dwrite(hdf5_dataset, hdf5_datatype, H5S_ALL, H5S_ALL,
                    H5P_DEFAULT, vel_buffer);
            H5Dclose(hdf5_dataset);
            H5Sclose(hdf5_dataspace);
            H5Tclose(hdf5_datatype);
            
            strcpy(buf, "Coordinates");
            rank=2;
            dims[1]=3;

            double* pos_buffer;
            pos_buffer = (double*) malloc ( dims[0] * dims[1] * sizeof(double) );

            for (int i = 0; i < dims[0]; i++ )
                for (int j = 0; j < dims[1]; j++ )
                {
                    int ipart = dims[1] * i + j;
                    pos_buffer[ipart] = (double) P[i].Pos[j];
                }
            
            hdf5_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
            hdf5_dataspace = H5Screate_simple(rank, dims, NULL);
            hdf5_dataset = H5Dcreate(hdf5_grp[type], buf, hdf5_datatype,
                    hdf5_dataspace, H5P_DEFAULT);
            H5Dwrite(hdf5_dataset, hdf5_datatype, H5S_ALL, H5S_ALL,
                    H5P_DEFAULT, pos_buffer);
            H5Dclose(hdf5_dataset);
            H5Sclose(hdf5_dataspace);
            H5Tclose(hdf5_datatype);

            strcpy(buf, "ParticleIDs");
            dims[1] = 1;
            rank=1;
            
            int32_t* ID_buffer;
            ID_buffer = (int32_t*) malloc ( dims[0] * dims[1] * sizeof(int32_t) );

            for (int i = 0; i < dims[0]; i++ )
                ID_buffer[i] = (int32_t) P[i].ID;

            hdf5_datatype = H5Tcopy(H5T_NATIVE_UINT32);
            hdf5_dataspace = H5Screate_simple(rank, dims, NULL);
            hdf5_dataset = H5Dcreate(hdf5_grp[type], buf, hdf5_datatype,
                    hdf5_dataspace, H5P_DEFAULT);
            H5Dwrite(hdf5_dataset, hdf5_datatype, H5S_ALL, H5S_ALL,
                    H5P_DEFAULT, ID_buffer);
            H5Dclose(hdf5_dataset);
            H5Sclose(hdf5_dataspace);
            H5Tclose(hdf5_datatype);

            strcpy(buf, "Type");
            dims[1] = 1;
            rank=1;
            
            float* type_buffer;
            type_buffer = (float*) malloc ( dims[0] * dims[1] * sizeof(float) );

            for (int i = 0; i < dims[0]; i++ )
                type_buffer[i] = (float) P[i].CMZ_type;

            hdf5_dataspace = H5Screate_simple(rank, dims, NULL);
            hdf5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);
            hdf5_dataset = H5Dcreate(hdf5_grp[type], buf, hdf5_datatype,
                    hdf5_dataspace, H5P_DEFAULT);
            H5Dwrite(hdf5_dataset, hdf5_datatype, H5S_ALL, H5S_ALL,
                    H5P_DEFAULT, type_buffer);
            H5Dclose(hdf5_dataset);
            H5Tclose(hdf5_datatype);
            H5Sclose(hdf5_dataspace);
        }
    }
    for(type = 5; type >= 0; type--)
        if(Header.npart[type] > 0)
            H5Gclose(hdf5_grp[type]);
    H5Gclose(hdf5_headergrp);
	H5Fclose(hdf5_file);
    return;
}

void write_header_attributes_in_hdf5(hid_t handle)
{
  hsize_t adim[1] = { 6 };
  hid_t hdf5_dataspace, hdf5_attribute;

  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
  hdf5_attribute = H5Acreate(handle, "NumPart_ThisFile", H5T_NATIVE_UINT64, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, Header.npart);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
  hdf5_attribute = H5Acreate(handle, "NumPart_Total", H5T_NATIVE_UINT64, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, Header.npartTotal);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
  hdf5_attribute = H5Acreate(handle, "MassTable", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, Header.mass);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "BoxSize", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &Header.BoxSize);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  double time = 0.0;
  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Time", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &time);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  double redshift = 0.0;
  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Redshift", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &redshift);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "NumFilesPerSnapshot", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &Header.num_files);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);
}

#endif
