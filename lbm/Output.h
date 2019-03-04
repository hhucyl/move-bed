#ifndef LBM_OUTPUT_H
#define LBM_OUTPUT_H

void Domain::WriteXDMF(char const * FileKey)
{
    String fn(FileKey);
    fn.append(".h5");
    hid_t     file_id;
    file_id = H5Fcreate(fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    size_t  Nx = Ndim(0);
    size_t  Ny = Ndim(1);
    size_t  Nz = Ndim(2);
    size_t Step = 1;
    size_t Nl = 1;
    for (size_t j=0;j<Nl;j++)
    {
        // Creating data sets
        double * Density   = new double[  Nx*Ny*Nz];
        double * Ga     = new double[  Nx*Ny*Nz];
        double * Overlap = new double [Nx*Ny*Nz];
        double * Vvec      = new double[3*Nx*Ny*Nz];
        double * BFvec      = new double[3*Nx*Ny*Nz];
        double * Vvecp      = new double[3*Nx*Ny*Nz];
        double * Vflbm      = new double[3*Nx*Ny*Nz];
        double * Ff = NULL;
        double * Fft = NULL;
        double *qq = NULL;
        if(Isq)  qq = new double[Nneigh*Nx*Ny*Nz];
        if(IsF)      Ff   = new double[Nneigh*Nx*Ny*Nz];
        if(IsF)      Fft   = new double[Nneigh*Nx*Ny*Nz];
        size_t i=0;
        for (size_t m=0;m<Nz;m+=Step)
        for (size_t l=0;l<Ny;l+=Step)
        for (size_t n=0;n<Nx;n+=Step)
        {
            double rho    = 0.0;
            double gamma  = 0.0;
            Vec3_t vel    = Vec3_t(0.0,0.0,0.0);
            Vec3_t velp    = Vec3_t(0.0,0.0,0.0);
            Vec3_t flbm    = Vec3_t(0.0,0.0,0.0);
            Vec3_t BF     = Vec3_t(0.0,0.0,0.0);
            double temp = 0.0;
            for (size_t ni=0;ni<Step;ni++)
            for (size_t li=0;li<Step;li++)
            for (size_t mi=0;mi<Step;mi++)
            {
                rho    += Rho    [n+ni][l+li][m+mi];
                temp    = IsSolid[n+ni][l+li][m+mi] ? 2.0: 0.0;
                gamma  += std::max(Gamma[n+ni][l+li][m+mi],temp);
                vel    += Vel    [n+ni][l+li][m+mi];
                BF    += BForce    [n+ni][l+li][m+mi];
                velp    += VelP    [n+ni][l+li][m+mi];
                flbm    += Flbm    [n+ni][l+li][m+mi];
            }
            rho  /= Step*Step*Step;
            gamma/= Step*Step*Step;
            vel  /= Step*Step*Step;
            velp  /= Step*Step*Step;
            flbm  /= Step*Step*Step;
            BF   /= Step*Step*Step;
            Ga   [i]  = (double) gamma;
            Density [i]  = (double) rho;            
            Vvec[3*i  ]  = (double) vel(0)*(1.0-Ga[i]);
            Vvec[3*i+1]  = (double) vel(1)*(1.0-Ga[i]);
            Vvec[3*i+2]  = (double) vel(2)*(1.0-Ga[i]);
            Vvecp[3*i  ]  = (double) velp(0);
            Vvecp[3*i+1]  = (double) velp(1);
            Vvecp[3*i+2]  = (double) velp(2);
            Vflbm[3*i  ]  = (double) flbm(0);
            Vflbm[3*i+1]  = (double) flbm(1);
            Vflbm[3*i+2]  = (double) flbm(2);
            BFvec[3*i ]   = (double) BF(0);
            BFvec[3*i+1]   = (double) BF(1);
            BFvec[3*i+2]   = (double) BF(2);
            if(IsF) 
            {
                for (size_t k=0; k<Nneigh; k++)
                {
                    Ff[Nneigh*i + k] = (double) F[n][l][m][k];
                }
            }
            if(IsFt) 
            {
                for (size_t k=0; k<Nneigh; k++)
                {
                    Fft[Nneigh*i + k] = (double) Ftemp[n][l][m][k];
                }
            }
            if(Isq) 
            {
                for (size_t k=0; k<Nneigh; k++)
                {
                    qq[Nneigh*i + k] = (double) q[n][l][m][k];
                }
            }
            i++;
        }
        double *Pposition = NULL;
        double *Ppositionb = NULL;
        double *PVeloc = NULL;
        double *PForce = NULL;
        double *PForceh = NULL;
        double *PW = NULL;
        double *PWb = NULL;
        int *Ptag = NULL;
        int *Plist = NULL;
        // int *PlistPP = NULL;
        // int *PlistPG = NULL;
        int NP = Particles.size()*2;
        int NL = ListofContacts.size()*2;
        // int NL = ListofContactsPP.size()*2;
        // int NLG = ListofContactsPG.size()*2;
        if(Particles.size()>0)
        {
            Pposition = new double[3*NP];
            Ppositionb = new double[3*NP];
            PVeloc = new double[3*NP];
            PForce = new double[3*NP];
            PForceh = new double[3*NP];
            PW = new double[3*NP];
            PWb = new double[3*NP];
            Ptag = new int[NP];
            Plist = new int[NL];
            // PlistPP = new int[NL];
            // PlistPG = new int[NLG];
            for(size_t ip=0; ip<Particles.size();++ip)
            {
                DEM::Disk *Pa = &Particles[ip];
                Ptag[ip] = 1;
                Pposition[3*ip] = Pa->X(0);
                Pposition[3*ip+1] = Pa->X(1);
                Pposition[3*ip+2] = Pa->X(2);
                Ppositionb[3*ip] = Pa->Xb(0);
                Ppositionb[3*ip+1] = Pa->Xb(1);
                Ppositionb[3*ip+2] = Pa->Xb(2);
                PVeloc[3*ip] = Pa->V(0);
                PVeloc[3*ip+1] = Pa->V(1);
                PVeloc[3*ip+2] = Pa->V(2);
                PW[3*ip] = Pa->W(0);
                PW[3*ip+1] = Pa->W(1);
                PW[3*ip+2] = Pa->W(2);
                PWb[3*ip] = Pa->Wb(0);
                PWb[3*ip+1] = Pa->Wb(1);
                PWb[3*ip+2] = Pa->Wb(2);
                PForce[3*ip] = Pa->Fc(0);
                PForce[3*ip+1] = Pa->Fc(1);
                PForce[3*ip+2] = Pa->Fc(2);
                PForceh[3*ip] = Pa->Fh(0);
                PForceh[3*ip+1] = Pa->Fh(1);
                PForceh[3*ip+2] = Pa->Fh(2);
            }
            int GN = Particles.size();
            for(int ip=0; ip<(int)GhostParticles.size();++ip)
            {
                DEM::Disk *Pa = &GhostParticles[ip];
                int ipp = ip+GN;
                Ptag[ipp] = Pa->Ghost ? 1.0 : -1.0;
                Pposition[3*ipp] = Pa->X(0);
                Pposition[3*ipp+1] = Pa->X(1);
                Pposition[3*ipp+2] = Pa->X(2);
                Ppositionb[3*ipp] = Pa->Xb(0);
                Ppositionb[3*ipp+1] = Pa->Xb(1);
                Ppositionb[3*ipp+2] = Pa->Xb(2);
                PVeloc[3*ipp] = Pa->V(0);
                PVeloc[3*ipp+1] = Pa->V(1);
                PVeloc[3*ipp+2] = Pa->V(2);
                PW[3*ipp] = Pa->W(0);
                PW[3*ipp+1] = Pa->W(1);
                PW[3*ipp+2] = Pa->W(2);
                PWb[3*ipp] = Pa->Wb(0);
                PWb[3*ipp+1] = Pa->Wb(1);
                PWb[3*ipp+2] = Pa->Wb(2);
                PForce[3*ipp] = Pa->Fc(0);
                PForce[3*ipp+1] = Pa->Fc(1);
                PForce[3*ipp+2] = Pa->Fc(2);
                PForceh[3*ipp] = Pa->Fh(0);
                PForceh[3*ipp+1] = Pa->Fh(1);
                PForceh[3*ipp+2] = Pa->Fh(2);
            }
            for(size_t il=0; il<ListofContacts.size(); ++il)
            {
                Plist[2*il] = ListofContacts[il].first;
                Plist[2*il+1] = ListofContacts[il].second;
            }
            // for(size_t il=0; il<ListofContactsPP.size(); ++il)
            // {
            //     PlistPP[2*il] = ListofContactsPP[il].first;
            //     PlistPP[2*il+1] = ListofContactsPP[il].second;
            // }
            // for(size_t il=0; il<ListofContactsPG.size(); ++il)
            // {
            //     PlistPG[2*il] = ListofContactsPG[il].first;
            //     PlistPG[2*il+1] = ListofContactsPG[il].second;
            // }
        }

        double *RWPpos = NULL;
        int Nrwp = RWParticles.size();
        if(RWParticles.size()>0)
        {
            RWPpos = new double[3*Nrwp];
            for(size_t i=0; i<RWParticles.size(); ++i)
            {
                RW::Particle *RWP = &RWParticles[i];
                RWPpos[3*i] = RWP->X(0);
                RWPpos[3*i+1] = RWP->X(1);
                RWPpos[3*i+2] = RWP->X(2);
            }
        }

        
        //Writing data to h5 file
        hsize_t dims[1];
        dims[0] = Nx*Ny*Nz;
        String dsname;
        dsname.Printf("Density_%d",j);
        H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,Density );
        if (j==0)
        {
            dsname.Printf("Gamma");
            H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,Ga   );
        }
        dims[0] = 3*Nx*Ny*Nz;
        dsname.Printf("Velocity_%d",j);
        H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,Vvec    );
        dims[0] = 3*Nx*Ny*Nz;
        dsname.Printf("Velocity_P_%d",j);
        H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,Vvecp    );
        dims[0] = 3*Nx*Ny*Nz;
        dsname.Printf("Flbm_%d",j);
        H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,Vflbm    );
        dims[0] = 3*Nx*Ny*Nz;
        dsname.Printf("BForce_%d",j);
        H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,BFvec    );
        dims[0] = Nneigh*Nx*Ny*Nz;
        
        if(IsF)
        {
            dsname.Printf("F_%d",j);
            H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,Ff    );
        }
        if(IsFt)
        {
            dsname.Printf("Ftemp_%d",j);
            H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,Fft    );
        }
        if(Isq)
        {
            dsname.Printf("q_%d",j);
            H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,qq    );
        }
        dims[0] = 1;
        int N[1];
        N[0] = Nx;
        dsname.Printf("Nx");
        H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,N);
        dims[0] = 1;
        N[0] = Ny;
        dsname.Printf("Ny");
        H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,N);
        dims[0] = 1;
        N[0] = Nz;
        dsname.Printf("Nz");
        H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,N);
        if(Particles.size()>0)
        {
            dims[0] = NP;
            dsname.Printf("PTag");        
            H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,Ptag);
            dims[0] = 3*NP;
            dsname.Printf("Pposition");        
            H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,Pposition);
            dsname.Printf("Ppositionb");        
            H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,Ppositionb);
            dsname.Printf("PVeloc");        
            H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,PVeloc);
            dsname.Printf("PW");        
            H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,PW);
            dsname.Printf("PWb");        
            H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,PWb);    
            dsname.Printf("PForce");        
            H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,PForce);
            dsname.Printf("PForceh");        
            H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,PForceh);  
            dims[0] = NL;
            dsname.Printf("PList");        
            H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,Plist);
            // dsname.Printf("PListPP");        
            // H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,PlistPP);
            // dims[0] = NLG;
            // dsname.Printf("PListPG");        
            // H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,PlistPG);    
           
        }
        if(RWParticles.size())
        {
            dims[0] = 3*Nrwp;
            dsname.Printf("RWPposition");        
            H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,RWPpos);
        }
        
        delete [] Density ;
        delete [] Ga   ;
        delete [] Vvec    ;
        delete [] Vvecp    ;
        delete [] Vflbm    ;
        delete [] BFvec  ;
        if(IsF) delete [] Ff; 
        if(IsFt) delete [] Fft; 
        if(Isq) delete [] qq; 
        delete [] Overlap;
        if(Particles.size()>0)
        {
            delete [] Ptag;
            delete [] Pposition;
            delete [] Ppositionb;
            delete [] PVeloc;
            delete [] PForce;
            delete [] PForceh;
            delete [] PW;
            delete [] PWb;
            // delete [] PlistPP;
            // delete [] PlistPG;
            delete [] Plist;
        }
        if(RWParticles.size()>0)
        {
            delete [] RWPpos;
        }
    }


    //Closing the file
    H5Fflush(file_id,H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);

    // Writing xmf file
    std::ostringstream oss;


    if (Nz==1)
    {
        oss << "<?xml version=\"1.0\" ?>\n";
        oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
        oss << "<Xdmf Version=\"2.0\">\n";
        oss << " <Domain>\n";
        oss << "   <Grid Name=\"mesh1\" GridType=\"Uniform\">\n";
        oss << "     <Topology TopologyType=\"2DCoRectMesh\" Dimensions=\"" << Ny << " " << Nx << "\"/>\n";
        oss << "     <Geometry GeometryType=\"ORIGIN_DXDY\">\n";
        oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"2\"> 0.0 0.0\n";
        oss << "       </DataItem>\n";
        oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"2\"> 1.0 1.0\n";
        oss << "       </DataItem>\n";
        oss << "     </Geometry>\n";
        for (size_t j=0;j<Nl;j++)
        {
        oss << "     <Attribute Name=\"Density_" << j << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Nx << " " << Ny << " " << Nz << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Density_" << j << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "     <Attribute Name=\"Velocity_" << j << "\" AttributeType=\"Vector\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Nx << " " << Ny << " " << Nz << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Velocity_" << j << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "     <Attribute Name=\"Velocity_P_" << j << "\" AttributeType=\"Vector\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Nx << " " << Ny << " " << Nz << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Velocity_P_" << j << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "     <Attribute Name=\"Flbm_" << j << "\" AttributeType=\"Vector\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Nx << " " << Ny << " " << Nz << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Flbm_" << j << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "     <Attribute Name=\"BForce_" << j << "\" AttributeType=\"Vector\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Nx << " " << Ny << " " << Nz << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/BForce_" << j << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        }
        oss << "     <Attribute Name=\"Gamma\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Nx << " " << Ny << " " << Nz << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Gamma\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "     <Attribute Name=\"Overlap\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Nx << " " << Ny << " " << Nz << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Overlap\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "   </Grid>\n";
        // oss << "   <Grid Name=\"DEM\" GridType=\"Uniform\">\n";
        // oss << "     <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << Particles.size() << "\"/>\n";
        // oss << "     <Geometry GeometryType=\"XYZ\">\n";
        // oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << Particles.size() << " 3\" >\n";
        // oss << "        " << fn.CStr() <<":/Pposition \n";
        // oss << "       </DataItem>\n";
        // oss << "     </Geometry>\n";
        // oss << "     <Attribute Name=\"PTag\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        // oss << "       <DataItem Dimensions=\"" << Particles.size() << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        // oss << "        " << fn.CStr() <<":/PTag \n";
        // oss << "       </DataItem>\n";
        // oss << "     </Attribute>\n";
        // oss << "     <Attribute Name=\"PForce\" AttributeType=\"Vector\" Center=\"Node\">\n";
        // oss << "       <DataItem Dimensions=\"" << Particles.size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        // oss << "        " << fn.CStr() <<":/PForce\n";
        // oss << "       </DataItem>\n";
        // oss << "     </Attribute>\n";
        // oss << "   </Grid>\n";
        oss << " </Domain>\n";
        oss << "</Xdmf>\n";
    }
    else
    {
        oss << "<?xml version=\"1.0\" ?>\n";
        oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
        oss << "<Xdmf Version=\"2.0\">\n";
        oss << " <Domain>\n";
        oss << "   <Grid Name=\"LBM_Mesh\" GridType=\"Uniform\">\n";
        oss << "     <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"" << Nz << " " << Ny << " " << Nx << "\"/>\n";
        oss << "     <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n";
        oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"3\"> 0.0 0.0 0.0\n";
        oss << "       </DataItem>\n";
        oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"3\"> " << Step*dx << " " << Step*dx  << " " << Step*dx  << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Geometry>\n";
        for (size_t j=0;j<Nl;j++)
        {
        oss << "     <Attribute Name=\"Density_" << j << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Nz << " " << Ny << " " << Nx << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Density_" << j << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "     <Attribute Name=\"Velocity_" << j << "\" AttributeType=\"Vector\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Nz << " " << Ny << " " << Nx << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Velocity_" << j << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "     <Attribute Name=\"Velocity_P_" << j << "\" AttributeType=\"Vector\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Nx << " " << Ny << " " << Nz << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Velocity_P_" << j << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "     <Attribute Name=\"BForce_" << j << "\" AttributeType=\"Vector\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Nz << " " << Ny << " " << Nx << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/BForce_" << j << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "     <Attribute Name=\"Flbm_" << j << "\" AttributeType=\"Vector\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Nx << " " << Ny << " " << Nz << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Flbm_" << j << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        
        }
        oss << "     <Attribute Name=\"Gamma\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Nz << " " << Ny << " " << Nx << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Gamma\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "   </Grid>\n";
        oss << " </Domain>\n";
        oss << "</Xdmf>\n";
    }
    fn = FileKey;
    fn.append(".xmf");
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close();

}

#endif

