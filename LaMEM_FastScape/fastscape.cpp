#include <iostream>
#include <fstream>
#include <random>

#include <stdio.h>
#include <string.h>

// LaMEM header file
#include "LaMEM.h"
#include "surf.h"
#include "scaling.h"
#include "JacRes.h"
#include "tssolve.h"
#include "advect.h"
#include "interpolate.h"
#include "fastscape.h"
#include "paraViewOutSurf.h"
#include "paraViewOutBin.h"

//---------------------------------------------------------------------------
PetscErrorCode FastScapeCreate(FastScapeLib *FSLib, FB *fb)
{
    PetscErrorCode ierr;
    PetscInt maxPhaseID;
    Scaling  *scal;
    FreeSurf       *surf;

    // access context
    surf    = FSLib->surf;
    scal    = FSLib->scal;

    //=================================================================================
    // Load data from .dat file
    //=================================================================================
    // initialize
    // non uniform grid
    FSLib->non_uniform_grid  =     0;
    // 2D grid
    FSLib->fs2D              =     0;
    // extend range & nodes
    FSLib->extendedNodes     =   101;
    FSLib->extendedRange     =   100 * scal->length_fs;
    // refine times & load refined grid
    FSLib->refine            =     1; 
    // max timestep
    FSLib->Max_dt            =  0.01 * scal->time_fs;
    // random noise
    FSLib->random_noise      =     1; 
    // sedimentation
    FSLib->setMarine         =     0; 
    // output information
    FSLib->surf_out_nstep    =     1;
    FSLib->vec_times         =     1;
    // total phases
    maxPhaseID = FSLib->surf->jr->dbm->numPhases-1;
 
    // load information from .dat file
    // setup block access mode                                                                                                                                      
    ierr = FBFindBlocks(fb, _REQUIRED_, "<FastScapeStart>", "<FastScapeEnd>");      CHKERRQ(ierr);

    if(fb->nblocks)
    {
        //-------------------------------
        // Grid information
        //-------------------------------
        // non uniform grid
        ierr         = getIntParam   (fb, _OPTIONAL_, "non_uniform_grid", &FSLib->non_uniform_grid, 1,        1);      CHKERRQ(ierr); // flag 
        // 2D grid
        ierr         = getIntParam   (fb, _OPTIONAL_, "fs2D",             &FSLib->fs2D,             1,        1);      CHKERRQ(ierr); // flag
        if( 1 == FSLib->fs2D )
        {
            ierr     = getScalarParam(fb, _REQUIRED_, "extendedRange",    &FSLib->extendedRange,    1, 1 / scal->length_fs);      CHKERRQ(ierr); // km (LaMEM) -> m (FastScape)
            
            if ( 0 == FSLib->non_uniform_grid)
            {
                ierr = getIntParam   (fb, _REQUIRED_, "extendedNodes",    &FSLib->extendedNodes,    1, 10000000);      CHKERRQ(ierr); // non-dimensional
            }
        }    
        // refined grid
        ierr         = getIntParam   (fb, _OPTIONAL_, "fs_refine",        &FSLib->refine,           1,     100);       CHKERRQ(ierr); // non-dimensional
       
        //===============================                                                                                                                                            
        // FastScape PARAMETER                                                                                                               
        //===============================
        // dt & boundary condition
        ierr = getScalarParam(fb, _REQUIRED_, "Max_dt",            &FSLib->Max_dt,           1, 1 / scal->time_fs);         CHKERRQ(ierr); // Myr (LaMEM) ->yr (FastScape)
        
        // bottom-right-top-left; 0 = reflective, 1 = fixed height boundary; When two reflective boundaris face each other they become cyclic
        ierr = getStringParam(fb, _REQUIRED_, "topo_boundary",     FSLib->FS_BC,                    "1111");         CHKERRQ(ierr); 
        
        // 1 -- boundary velocity == 0; 0 -- boundary velocity from LaMEM
        ierr = getStringParam(fb, _REQUIRED_, "vel_boundary",      FSLib->FS_VELBC,                 "1111");         CHKERRQ(ierr); 
        
        // random noise
        ierr = getIntParam   (fb, _REQUIRED_, "random_noise",      &FSLib->random_noise,     1,         1);         CHKERRQ(ierr); 

        // sedimentation phase
        ierr = getIntParam   (fb, _REQUIRED_, "sed_phases",        &FSLib->sedPhases,        1,   maxPhaseID);         CHKERRQ(ierr); // non-dimensional

        surf->phase = FSLib->sedPhases;
        //-------------------------------
        // Erosion process
        //-------------------------------
        // kf, kd can be set as an array
        ierr = getScalarParam(fb, _REQUIRED_, "kf",                &FSLib->kf,               1,          1.0);         CHKERRQ(ierr); // m/yr
        ierr = getScalarParam(fb, _REQUIRED_, "kfsed",             &FSLib->kfsed,            1,          1.0);         CHKERRQ(ierr); // m/yr
        ierr = getScalarParam(fb, _REQUIRED_, "m",                 &FSLib->m,                1,          1.0);         CHKERRQ(ierr); // non-dimensional
        ierr = getScalarParam(fb, _REQUIRED_, "n",                 &FSLib->n,                1,          1.0);         CHKERRQ(ierr); // non-dimensional
        ierr = getScalarParam(fb, _REQUIRED_, "kd",                &FSLib->kd,               1,          1.0);         CHKERRQ(ierr); // m/yr
        ierr = getScalarParam(fb, _REQUIRED_, "kdsed",             &FSLib->kdsed,            1,          1.0);         CHKERRQ(ierr); // m/yr
        ierr = getScalarParam(fb, _REQUIRED_, "g",                 &FSLib->g,                1,          1.0);         CHKERRQ(ierr); // non-dimensional
        ierr = getScalarParam(fb, _REQUIRED_, "gsed",              &FSLib->gsed,             1,          1.0);         CHKERRQ(ierr); // non-dimensional
        ierr = getScalarParam(fb, _REQUIRED_, "p",                 &FSLib->p,                1,          1.0);         CHKERRQ(ierr); // non-dimensional

        //-------------------------------
        // Sedimentation process
        //-------------------------------    
        ierr = getIntParam   (fb, _OPTIONAL_, "setMarine",         &FSLib->setMarine,        1,            1);         CHKERRQ(ierr); // flag
        if(1 == FSLib->setMarine)
        {
            ierr = getScalarParam(fb, _REQUIRED_, "sealevel",      &FSLib->sealevel,         1,  1 / scal->length_fs);          CHKERRQ(ierr); // m
            ierr = getScalarParam(fb, _REQUIRED_, "poroSilt",      &FSLib->poro_silt,        1,         1.0);          CHKERRQ(ierr); // non-dimensional
            ierr = getScalarParam(fb, _REQUIRED_, "poroSand",      &FSLib->poro_sand,        1,         1.0);          CHKERRQ(ierr); // non-dimensional
            ierr = getScalarParam(fb, _REQUIRED_, "zporoSilt",     &FSLib->zporo_silt,       1,         1.0);          CHKERRQ(ierr); 
            ierr = getScalarParam(fb, _REQUIRED_, "zporoSand",     &FSLib->zporo_sand,       1,         1.0);          CHKERRQ(ierr); 
            ierr = getScalarParam(fb, _REQUIRED_, "ratio",         &FSLib->ratio,            1,         1.0);          CHKERRQ(ierr); // non-dimensional
            ierr = getScalarParam(fb, _REQUIRED_, "Lsolve",        &FSLib->Lsolve,           1,         1.0);          CHKERRQ(ierr); // m
            ierr = getScalarParam(fb, _REQUIRED_, "kdsSilt",       &FSLib->kds_silt,         1,         1.0);          CHKERRQ(ierr); // m2/yr
            ierr = getScalarParam(fb, _REQUIRED_, "kdsSand",       &FSLib->kds_sand,         1,         1.0);          CHKERRQ(ierr); // m2/yr
        }

        //-------------------------------
        // Output information
        //-------------------------------
        ierr = getIntParam   (fb, _OPTIONAL_, "surf_out_nstep",    &FSLib->surf_out_nstep,   1,        1e6);           CHKERRQ(ierr); // non-dimensional
        ierr = getScalarParam(fb, _OPTIONAL_, "vec_times",         &FSLib->vec_times,        1,        1.0);           CHKERRQ(ierr); // non-dimensional
    }

    ierr = FBFreeBlocks(fb);       CHKERRQ(ierr);

    //=================================================================================
    // Load grid information from LaMEM
    //=================================================================================
    // load nx, ny, rangeX, rangeY, rangeZ
    ierr = FastScapeLoadGridInf(FSLib);  CHKERRQ(ierr);
    ierr = FastScapeCreateSurfaceGrid(FSLib, 1);  CHKERRQ(ierr);

    // save nx & ny
    if(1 == FSLib->fs2D)
    {
        if(1 == FSLib->refine)
        {
            FSLib->nx_solve    = FSLib->extendedXNodes;
            FSLib->ny_solve    = FSLib->extendedYNodes;
        }
        else
        {
            FSLib->nx_solve    = FSLib->etRefineXNodes;
            FSLib->ny_solve    = FSLib->etRefineYNodes;
        }
    }
    else
    {
        if(1 == FSLib->refine)
        {
            if( 0 == FSLib->non_uniform_grid)
            {
                FSLib->nx_solve    = FSLib->nx_fs;
                FSLib->ny_solve    = FSLib->ny_fs;
            }
            else
            {
                FSLib->nx_solve    = FSLib->msx_fs.nnodes_nug;
                FSLib->ny_solve    = FSLib->msy_fs.nnodes_nug;                
            }
        }
        else
        {
            FSLib->nx_solve    = FSLib->nx_refine;
            FSLib->ny_solve    = FSLib->ny_refine;
        }        
    }
    
    FSLib->nodes_solve = FSLib->nx_solve * FSLib->ny_solve;

    //=================================================================================
    // Visualization
    //=================================================================================
    PetscPrintf(PETSC_COMM_WORLD, "FastScape parameters: \n");
    // LaMEM grid
    PetscPrintf(PETSC_COMM_WORLD, "    Original grid: \n");  
    PetscPrintf(PETSC_COMM_WORLD, "    [nodeX, nodeY]        : [%d, %d]\n",     FSLib->nx_fs, FSLib->ny_fs);        
    // non uniform grid
    if( 1 == FSLib->non_uniform_grid )
    {
        PetscPrintf(PETSC_COMM_WORLD, "    Non unifrom grid: \n");  
        PetscPrintf(PETSC_COMM_WORLD, "    [nodeX, nodeY]        : [%d, %d]\n",     FSLib->msx_fs.nnodes_nug, FSLib->msy_fs.nnodes_nug);         
    }
    // 2D extended grid
    if( 1 == FSLib->fs2D )  
    {
        PetscPrintf(PETSC_COMM_WORLD, "    Extended  grid: \n");
        PetscPrintf(PETSC_COMM_WORLD, "    [rangeX,rangeY]       : [%g, %g] %s\n",   FSLib->extendedXRange / scal->length_fs, FSLib->extendedYRange / scal->length_fs, scal->lbl_length);
        PetscPrintf(PETSC_COMM_WORLD, "    [nodeX, nodeY]        : [%d, %d]\n",      FSLib->extendedXNodes, FSLib->extendedYNodes);
    }
    // refined grid
    if( 1 < FSLib->refine ) 
    {
        PetscPrintf(PETSC_COMM_WORLD, "    Refined grid: \n");
        PetscPrintf(PETSC_COMM_WORLD, "    Refined times         : %d\n",       FSLib->refine);   
        // 2D
        if( 1 == FSLib->fs2D)
        {        
            PetscPrintf(PETSC_COMM_WORLD, "    [nodeX, nodeY]        : [%d, %d]\n", FSLib->etRefineXNodes, FSLib->etRefineYNodes);        
        }
        // 3D
        else
        {
            PetscPrintf(PETSC_COMM_WORLD, "    [nodeX, nodeY]        : [%d, %d]\n", FSLib->nx_refine, FSLib->ny_refine);  
        }
    }

    // surface process parameter
    PetscPrintf(PETSC_COMM_WORLD, "  Surface process: \n");  
    PetscPrintf(PETSC_COMM_WORLD, "    Max timestep          : %g %s\n",        FSLib->Max_dt / scal->time_fs, scal->lbl_time);
    PetscPrintf(PETSC_COMM_WORLD, "    Topography boundary   : %s\n",           FSLib->FS_BC);
    PetscPrintf(PETSC_COMM_WORLD, "    Velocity boundary     : %s\n",           FSLib->FS_VELBC);
    PetscPrintf(PETSC_COMM_WORLD, "    Sedimentation phase   : %d\n",           FSLib->sedPhases);    
    PetscPrintf(PETSC_COMM_WORLD, "    SPL: \n");   
    PetscPrintf(PETSC_COMM_WORLD, "      Kf                  : %g\n",           FSLib->kf);
    PetscPrintf(PETSC_COMM_WORLD, "      Kfsed               : %g\n",           FSLib->kfsed);
    PetscPrintf(PETSC_COMM_WORLD, "      m                   : %g\n",           FSLib->m);    
    PetscPrintf(PETSC_COMM_WORLD, "      n                   : %g\n",           FSLib->n);
    PetscPrintf(PETSC_COMM_WORLD, "    Hillslope process: \n"); 
    PetscPrintf(PETSC_COMM_WORLD, "      Kd                  : %g\n",           FSLib->kd);
    PetscPrintf(PETSC_COMM_WORLD, "      Kdsed               : %g\n",           FSLib->kdsed);
    PetscPrintf(PETSC_COMM_WORLD, "      g                   : %g\n",           FSLib->g);    
    PetscPrintf(PETSC_COMM_WORLD, "      gsed                : %g\n",           FSLib->gsed);
    PetscPrintf(PETSC_COMM_WORLD, "      p                   : %g\n",           FSLib->p);    

    if(FSLib->setMarine == 1)
    {
        PetscPrintf(PETSC_COMM_WORLD, "    Marine process: \n");        
        PetscPrintf(PETSC_COMM_WORLD, "      sealevel            : %g %s\n",    FSLib->sealevel / scal->length_fs, scal->lbl_length);
        PetscPrintf(PETSC_COMM_WORLD, "      poro_silt           : %g\n",       FSLib->poro_silt);
        PetscPrintf(PETSC_COMM_WORLD, "      poro_sand           : %g\n",       FSLib->poro_sand);
        PetscPrintf(PETSC_COMM_WORLD, "      zporo_silt          : %g\n",       FSLib->zporo_silt);
        PetscPrintf(PETSC_COMM_WORLD, "      zporo_sand          : %g\n",       FSLib->zporo_sand);
        PetscPrintf(PETSC_COMM_WORLD, "      ratio               : %g\n",       FSLib->ratio);
        PetscPrintf(PETSC_COMM_WORLD, "      L                   : %g\n",       FSLib->Lsolve);
        PetscPrintf(PETSC_COMM_WORLD, "      kds_silt            : %g\n",       FSLib->kds_silt);
        PetscPrintf(PETSC_COMM_WORLD, "      kds_sand            : %g\n",       FSLib->kds_sand);        
    }
    PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");

    PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FastScapeCreateData(FastScapeLib *FSLib)
{
    FreeSurf       *surf;

    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    PetscInt nodes, nodes_nug, nodes_ori;

    // access context
    surf = FSLib->surf;

    // nodes
    nodes     = FSLib->nodes_solve;
    nodes_nug = FSLib->fsX.nodes * FSLib->fsY.nodes;
    nodes_ori = FSLib->nx_fs * FSLib->ny_fs;

    // non_uniform_grid
    ierr = VecCreate(PETSC_COMM_WORLD, &FSLib->gtopo_nug);            CHKERRQ(ierr);
    ierr = VecSetSizes(FSLib->gtopo_nug, nodes_nug, nodes_nug);               CHKERRQ(ierr);  
    ierr = VecSetFromOptions(FSLib->gtopo_nug);                       CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD, &FSLib->vx_nug);               CHKERRQ(ierr);
    ierr = VecSetSizes(FSLib->vx_nug, nodes_nug, nodes_nug);                  CHKERRQ(ierr);  
    ierr = VecSetFromOptions(FSLib->vx_nug);                          CHKERRQ(ierr);
    
    ierr = VecCreate(PETSC_COMM_WORLD, &FSLib->vy_nug);               CHKERRQ(ierr);
    ierr = VecSetSizes(FSLib->vy_nug, nodes_nug, nodes_nug);                  CHKERRQ(ierr);  
    ierr = VecSetFromOptions(FSLib->vy_nug);                          CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD, &FSLib->vz_nug);               CHKERRQ(ierr);
    ierr = VecSetSizes(FSLib->vz_nug, nodes_nug, nodes_nug);                  CHKERRQ(ierr);  
    ierr = VecSetFromOptions(FSLib->vz_nug);                          CHKERRQ(ierr);

    // vx & vy &vz
    ierr = DMCreateGlobalVector(surf->DA_SURF, &FSLib->vz_collect);   CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(surf->DA_SURF, &FSLib->vx_collect);   CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(surf->DA_SURF, &FSLib->vy_collect);   CHKERRQ(ierr);

    // topography 
    ierr = VecCreate(PETSC_COMM_WORLD, &FSLib->gtopo_fs);             CHKERRQ(ierr);
    ierr = VecSetSizes(FSLib->gtopo_fs, nodes_ori, nodes_ori);  CHKERRQ(ierr);  
    ierr = VecSetFromOptions(FSLib->gtopo_fs);                        CHKERRQ(ierr);

    // refined part
    ierr = VecCreate(PETSC_COMM_WORLD, &FSLib->gtopo_refine);         CHKERRQ(ierr);
    ierr = VecSetSizes(FSLib->gtopo_refine, nodes, nodes);            CHKERRQ(ierr);  
    ierr = VecSetFromOptions(FSLib->gtopo_refine);                    CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD, &FSLib->vx_refine);            CHKERRQ(ierr);
    ierr = VecSetSizes(FSLib->vx_refine, nodes, nodes);               CHKERRQ(ierr);  
    ierr = VecSetFromOptions(FSLib->vx_refine);  

    ierr = VecCreate(PETSC_COMM_WORLD, &FSLib->vy_refine);            CHKERRQ(ierr);
    ierr = VecSetSizes(FSLib->vy_refine, nodes, nodes);               CHKERRQ(ierr);  
    ierr = VecSetFromOptions(FSLib->vy_refine);  

    ierr = VecCreate(PETSC_COMM_WORLD, &FSLib->vz_refine);            CHKERRQ(ierr);
    ierr = VecSetSizes(FSLib->vz_refine, nodes, nodes);               CHKERRQ(ierr);  
    ierr = VecSetFromOptions(FSLib->vz_refine);                       CHKERRQ(ierr);

    // extended part
    ierr = VecCreate(PETSC_COMM_WORLD, &FSLib->gtopo_extend);         CHKERRQ(ierr);
    ierr = VecSetSizes(FSLib->gtopo_extend, nodes, nodes);            CHKERRQ(ierr);  
    ierr = VecSetFromOptions(FSLib->gtopo_extend);                    CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD, &FSLib->vx_extend);            CHKERRQ(ierr);
    ierr = VecSetSizes(FSLib->vx_extend, nodes, nodes);               CHKERRQ(ierr);  
    ierr = VecSetFromOptions(FSLib->vx_extend);                       CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD, &FSLib->vy_extend);            CHKERRQ(ierr);
    ierr = VecSetSizes(FSLib->vy_extend, nodes, nodes);               CHKERRQ(ierr);  
    ierr = VecSetFromOptions(FSLib->vy_extend);                       CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD, &FSLib->vz_extend);            CHKERRQ(ierr);
    ierr = VecSetSizes(FSLib->vz_extend, nodes, nodes);               CHKERRQ(ierr);  
    ierr = VecSetFromOptions(FSLib->vz_extend);                       CHKERRQ(ierr);

    // extended part after refinement
    ierr = VecCreate(PETSC_COMM_WORLD, &FSLib->gtopo_et_refine);      CHKERRQ(ierr);
    ierr = VecSetSizes(FSLib->gtopo_et_refine, nodes, nodes);         CHKERRQ(ierr);  
    ierr = VecSetFromOptions(FSLib->gtopo_et_refine);                 CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD, &FSLib->vx_et_refine);         CHKERRQ(ierr);
    ierr = VecSetSizes(FSLib->vx_et_refine, nodes, nodes);            CHKERRQ(ierr);  
    ierr = VecSetFromOptions(FSLib->vx_et_refine);                    CHKERRQ(ierr);
        
    ierr = VecCreate(PETSC_COMM_WORLD, &FSLib->vy_et_refine);         CHKERRQ(ierr);
    ierr = VecSetSizes(FSLib->vy_et_refine, nodes, nodes);            CHKERRQ(ierr);  
    ierr = VecSetFromOptions(FSLib->vy_et_refine);                    CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD, &FSLib->vz_et_refine);         CHKERRQ(ierr);
    ierr = VecSetSizes(FSLib->vz_et_refine, nodes, nodes);            CHKERRQ(ierr);  
    ierr = VecSetFromOptions(FSLib->vz_et_refine);                    CHKERRQ(ierr);
    
    // FastScape solution
    ierr = VecCreate(PETSC_COMM_WORLD, &FSLib->silt_fraction);        CHKERRQ(ierr);
    ierr = VecSetSizes(FSLib->silt_fraction, nodes, nodes);           CHKERRQ(ierr);  
    ierr = VecSetFromOptions(FSLib->silt_fraction);                   CHKERRQ(ierr);      

    ierr = VecCreate(PETSC_COMM_WORLD, &FSLib->basement);             CHKERRQ(ierr);
    ierr = VecSetSizes(FSLib->basement, nodes, nodes);                CHKERRQ(ierr);  
    ierr = VecSetFromOptions(FSLib->basement);                        CHKERRQ(ierr);  

    ierr = VecCreate(PETSC_COMM_WORLD, &FSLib->total_erosion);        CHKERRQ(ierr);
    ierr = VecSetSizes(FSLib->total_erosion, nodes, nodes);           CHKERRQ(ierr);  
    ierr = VecSetFromOptions(FSLib->total_erosion);                   CHKERRQ(ierr);  

    ierr = VecCreate(PETSC_COMM_WORLD, &FSLib->drainage_area);        CHKERRQ(ierr);
    ierr = VecSetSizes(FSLib->drainage_area, nodes, nodes);           CHKERRQ(ierr);  
    ierr = VecSetFromOptions(FSLib->drainage_area);                   CHKERRQ(ierr);  

    ierr = VecCreate(PETSC_COMM_WORLD, &FSLib->erosion_rate);         CHKERRQ(ierr);
    ierr = VecSetSizes(FSLib->erosion_rate, nodes, nodes);            CHKERRQ(ierr);  
    ierr = VecSetFromOptions(FSLib->erosion_rate);                    CHKERRQ(ierr);  

    ierr = VecCreate(PETSC_COMM_WORLD, &FSLib->slope);                CHKERRQ(ierr);
    ierr = VecSetSizes(FSLib->slope, nodes, nodes);                   CHKERRQ(ierr);  
    ierr = VecSetFromOptions(FSLib->slope);                           CHKERRQ(ierr);  

    ierr = VecCreate(PETSC_COMM_WORLD, &FSLib->curvature);            CHKERRQ(ierr);
    ierr = VecSetSizes(FSLib->curvature, nodes, nodes);               CHKERRQ(ierr);  
    ierr = VecSetFromOptions(FSLib->curvature);                       CHKERRQ(ierr);  

    ierr = VecCreate(PETSC_COMM_WORLD, &FSLib->chi);                  CHKERRQ(ierr);
    ierr = VecSetSizes(FSLib->chi, nodes, nodes);                     CHKERRQ(ierr);  
    ierr = VecSetFromOptions(FSLib->chi);                             CHKERRQ(ierr);  

    ierr = VecCreate(PETSC_COMM_WORLD, &FSLib->catchment);            CHKERRQ(ierr);
    ierr = VecSetSizes(FSLib->catchment, nodes, nodes);               CHKERRQ(ierr);  
    ierr = VecSetFromOptions(FSLib->catchment);                       CHKERRQ(ierr);  

    ierr = VecCreate(PETSC_COMM_WORLD, &FSLib->lake_depth);           CHKERRQ(ierr);
    ierr = VecSetSizes(FSLib->lake_depth, nodes, nodes);              CHKERRQ(ierr);  
    ierr = VecSetFromOptions(FSLib->lake_depth);                      CHKERRQ(ierr);      

    PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FastScapeLoadGridInf(FastScapeLib *FSLib)
{

    PetscErrorCode ierr;
    PetscScalar bx, by, bz, ex, ey, ez;

    // load global nx, ny, rangeX, rangeY, rangeZ
    FDSTAG   *fs;    
    Scaling  *scal;
    fs    = FSLib->surf->jr->fs;
    scal  = fs->scal;

    // range X, Y, Z   
    ierr  = FDSTAGGetGlobalBox(fs, &bx, &by, &bz, &ex, &ey, &ez); CHKERRQ(ierr);

    FSLib->rangeX_begin = bx * scal->length; // km (LaMEM) in GEO
    FSLib->rangeX_end   = ex * scal->length;
    FSLib->rangeY_begin = by * scal->length;
    FSLib->rangeY_end   = ey * scal->length;
    FSLib->rangeZ_begin = bz * scal->length;   
    FSLib->rangeZ_end   = ez * scal->length;

    FSLib->rangeX       = (FSLib->rangeX_end - FSLib->rangeX_begin) * scal->length_fs; //(km) in LaMEM to (m) in FastScape (GEO)
    FSLib->rangeY       = (FSLib->rangeY_end - FSLib->rangeY_begin) * scal->length_fs;

    // original nx, ny
    FSLib->nx_fs        = fs->dsx.tnods;
    FSLib->ny_fs        = fs->dsy.tnods;
    FSLib->fsX.nodes    = FSLib->nx_fs;
    FSLib->fsY.nodes    = FSLib->ny_fs;

    // create new nodes & range
    if(0 == FSLib->non_uniform_grid)
    {
        // 3D refine
        if(1 < FSLib->refine && 0 == FSLib->fs2D)
        {
            // refined grid in FastScaoe
            FSLib->nx_refine = (FSLib->nx_fs - 1) * FSLib->refine + 1;
            FSLib->ny_refine = (FSLib->ny_fs - 1) * FSLib->refine + 1;
            FSLib->fsX.nodes_refine      = FSLib->nx_refine;
            FSLib->fsY.nodes_refine      = FSLib->ny_refine;
        }
        // 2D
        if(1 == FSLib->fs2D)
        {
            // 2D No Refine
            // extend grid in FastScape
            // extend in rangeX
            if(FSLib->rangeX > FSLib->rangeY)
            {
                FSLib->extendedXRange = FSLib->rangeX;
                FSLib->extendedYRange = FSLib->extendedRange;
                FSLib->extendedXNodes = FSLib->nx_fs;
                FSLib->extendedYNodes = FSLib->extendedNodes;
                FSLib->extendedX      = 0;
                FSLib->extendedY      = 1;
            }
            // extend in rangeY    
            else
            {
                FSLib->extendedXRange = FSLib->extendedRange;
                FSLib->extendedYRange = FSLib->rangeY;
                FSLib->extendedXNodes = FSLib->extendedNodes;
                FSLib->extendedYNodes = FSLib->ny_fs;
                FSLib->extendedX      = 1;
                FSLib->extendedY      = 0;
            }
 
            FSLib->fsX.nodes_extend      = FSLib->extendedXNodes;
            FSLib->fsY.nodes_extend      = FSLib->extendedYNodes;

            // extended grid after refinement in FastScape
            if(1 < FSLib->refine)
            {
                FSLib->etRefineXNodes = (FSLib->extendedXNodes - 1) * FSLib->refine + 1;
                FSLib->etRefineYNodes = (FSLib->extendedYNodes - 1) * FSLib->refine + 1;

                FSLib->fsX.nodes_refine      = FSLib->etRefineXNodes;
                FSLib->fsY.nodes_refine      = FSLib->etRefineYNodes;
            }
        }
    }
    else
    {
        // original grid
        // setting bias-flag, minimum grid spacing, grid nodes
        // x-direction   
        ierr = FSLoadNonUniformGrid(&FSLib->msx_fs, FSLib->rangeX_end / scal->length, fs->scal);  CHKERRQ(ierr);
        // y-direction
        ierr = FSLoadNonUniformGrid(&FSLib->msy_fs, FSLib->rangeY_end / scal->length, fs->scal);  CHKERRQ(ierr);
        FSLib->fsX.nodes      = FSLib->msx_fs.nnodes_nug;
        FSLib->fsY.nodes      = FSLib->msy_fs.nnodes_nug;

        // 3D refine
        if(1 < FSLib->refine && 0 == FSLib->fs2D)
        {
            // refined grid in FastScaoe
            FSLib->nx_refine = (FSLib->fsX.nodes - 1) * FSLib->refine + 1;
            FSLib->ny_refine = (FSLib->fsY.nodes - 1) * FSLib->refine + 1;
            FSLib->fsX.nodes_refine      = FSLib->nx_refine;
            FSLib->fsY.nodes_refine      = FSLib->ny_refine;
        }
        // 2D
        if(1 == FSLib->fs2D)
        {
            // 2D No Refine
            // extend grid in FastScape
            // extend in rangeX
            if(FSLib->rangeX > FSLib->rangeY)
            {
                FSLib->extendedXRange = FSLib->rangeX;
                FSLib->extendedYRange = FSLib->extendedRange;
                FSLib->extendedXNodes = FSLib->fsX.nodes;
                FSLib->extendedYNodes = (PetscInt)(FSLib->extendedYRange / scal->length_fs / FSLib->msx_fs.min_spacing) + 2;           
                FSLib->extendedX      = 0;
                FSLib->extendedY      = 1;
            }
            // extend in rangeY    
            else
            {
                FSLib->extendedXRange = FSLib->extendedRange;
                FSLib->extendedYRange = FSLib->rangeY;
                FSLib->extendedXNodes = (PetscInt)(FSLib->extendedXRange / scal->length_fs / FSLib->msy_fs.min_spacing) + 2;
                FSLib->extendedYNodes = FSLib->fsY.nodes;
                FSLib->extendedX      = 1;
                FSLib->extendedY      = 0;
            }
            FSLib->fsX.nodes_extend      = FSLib->extendedXNodes;
            FSLib->fsY.nodes_extend      = FSLib->extendedYNodes;

            // extended grid after refinement in FastScape
            if(1 < FSLib->refine)
            {
                FSLib->etRefineXNodes = (FSLib->fsX.nodes_extend - 1) * FSLib->refine + 1;
                FSLib->etRefineYNodes = (FSLib->fsY.nodes_extend - 1) * FSLib->refine + 1;
                FSLib->fsX.nodes_refine      = FSLib->etRefineXNodes;
                FSLib->fsY.nodes_refine      = FSLib->etRefineYNodes;
            }
        }
    }

    PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FSLoadNonUniformGrid(MeshSeg1DFS *ms_fs, PetscScalar xend, Scaling *scal)
{
    PetscInt i;
    PetscScalar begSz, endSz, avgSz, bias;

    // initialize 
    ms_fs->bias = 0;
    ms_fs->xstart[ms_fs->nsegs] = xend;

    // setting bias-flag, minimum grid spacing, grid nodes
    // bias == 1.0 ? 
    for(i = 0; i < ms_fs->nsegs; i++)
    {
        if(ms_fs->biases[i] != 1.0)
        {
            ms_fs->bias = 1; break;
        }
    }

    // obtain the minimum grid spacing    
    for(i = 0; i < ms_fs->nsegs; i++)
    {
        if(0 == ms_fs->bias)
        {
            ms_fs->grid_spacing_min[i]  = (ms_fs->xstart[i+1] - ms_fs->xstart[i]) * scal->length / (ms_fs->istart[i+1] - ms_fs->istart[i]);
            ms_fs->grid_spacing_max[i]  = ms_fs->grid_spacing_min[i];
        }
        else
        {      
       
            // bias (last to first cell size ratio > 1 -> growing)
            bias  = ms_fs->biases[i];
        
            // average cell size
            avgSz = (ms_fs->xstart[i+1] - ms_fs->xstart[i]) * scal->length / (ms_fs->istart[i+1] - ms_fs->istart[i]);

            // cell size limits
            begSz = 2.0*avgSz/(1.0 + bias);
            endSz = bias*begSz;

            ms_fs->grid_spacing_min[i] = (begSz < endSz)? begSz : endSz;   
            ms_fs->grid_spacing_max[i] = (begSz > endSz)? begSz : endSz;                 
        }
    }    
    ms_fs->min_spacing = *min_element(ms_fs->grid_spacing_min, ms_fs->grid_spacing_min + ms_fs->nsegs);
    ms_fs->max_spacing = *max_element(ms_fs->grid_spacing_max, ms_fs->grid_spacing_max + ms_fs->nsegs);

    // get new spacing & nodes
    ms_fs->nnodes_nug  = floor( (ms_fs->xstart[ms_fs->nsegs] - ms_fs->xstart[0]) * scal->length / ms_fs->min_spacing ) + 2; 

    PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FastScapeCreateSurfaceGrid(FastScapeLib *FSLib, PetscInt mode)
{
    //mode = 1, normal run (create), = 2, restart, = 0, update range
    FSGrid   *fsX;
    FSGrid   *fsY;
    Scaling  *scal;
    TSSol    *ts;
    JacRes   *jr;
    PetscInt i, j, step_fs;
    PetscErrorCode ierr;

    fsX     = &FSLib->fsX;
    fsY     = &FSLib->fsY;
    scal    = FSLib->scal;
    jr      = FSLib->jr;
    ts      = jr->ts;
    step_fs = ts->istep;

    // original fastscape grid
    if(2 == mode || 1 == mode)
    {
        fsX->ncoor = (PetscScalar *) malloc(fsX->nodes * sizeof(PetscScalar));   
        fsY->ncoor = (PetscScalar *) malloc(fsY->nodes * sizeof(PetscScalar));   
    }

    fsX->dx = FSLib->rangeX / scal->length_fs / (fsX->nodes - 1);
    fsY->dx = FSLib->rangeY / scal->length_fs / (fsY->nodes - 1);
    
    // x-direction
    for(i = 0; i < fsX->nodes; i++)
    {
        fsX->ncoor[i] = FSLib->rangeX_begin + fsX->dx * i;
    }
    // y-direction
    for(j = 0; j < fsY->nodes; j++)
    {
        fsY->ncoor[j] = FSLib->rangeY_begin + fsY->dx * j;
    }

    if(0 == FSLib->fs2D)
    {
        // 3D refined grid
        if( 1 < FSLib->refine)
        {
            if(2 == mode || 0 == step_fs)
            {
                fsX->ncoor_refine = (PetscScalar *) malloc(fsX->nodes_refine * sizeof(PetscScalar));   
                fsY->ncoor_refine = (PetscScalar *) malloc(fsY->nodes_refine * sizeof(PetscScalar));          
            }

            fsX->dx_refine = FSLib->rangeX / scal->length_fs  / (fsX->nodes_refine - 1);
            fsY->dx_refine = FSLib->rangeY / scal->length_fs  / (fsY->nodes_refine - 1);

            // x-direction
            for(i = 0; i < fsX->nodes_refine; i++)
            {
                fsX->ncoor_refine[i] = FSLib->rangeX_begin + fsX->dx_refine * i;
            }
            // y-direction
            for(j = 0; j < fsY->nodes_refine; j++)
            {
                fsY->ncoor_refine[j] = FSLib->rangeY_begin + fsY->dx_refine * j;
            }        
        }
    }
    else
    {
        // 2D grid 
        if(2 == mode || 0 == step_fs)
        {
            fsX->ncoor_extend = (PetscScalar *) malloc(fsX->nodes_extend * sizeof(PetscScalar));   
            fsY->ncoor_extend = (PetscScalar *) malloc(fsY->nodes_extend * sizeof(PetscScalar)); 
        }

        fsX->dx_extend = FSLib->extendedXRange / scal->length_fs  / (fsX->nodes_extend - 1);
        fsY->dx_extend = FSLib->extendedYRange / scal->length_fs  / (fsY->nodes_extend - 1); 

        for(i = 0; i < fsX->nodes_extend; i++)
        {
            fsX->ncoor_extend[i] = FSLib->rangeX_begin + fsX->dx_extend * i;
        }
        // y-direction
        for(j = 0; j < fsY->nodes_extend; j++)
        {
            fsY->ncoor_extend[j] = FSLib->rangeY_begin + fsY->dx_extend * j;
        }        

        // 2D refined grid 
        if(1 < FSLib->refine)
        {
            if(2 == mode || 0 == step_fs)
            {
                fsX->ncoor_refine = (PetscScalar *) malloc(fsX->nodes_refine * sizeof(PetscScalar));   
                fsY->ncoor_refine = (PetscScalar *) malloc(fsY->nodes_refine * sizeof(PetscScalar));       
            }

            fsX->dx_refine = FSLib->extendedXRange / scal->length_fs  / (fsX->nodes_refine - 1);
            fsY->dx_refine = FSLib->extendedYRange / scal->length_fs  / (fsY->nodes_refine - 1);

            // x-direction
            for(i = 0; i < fsX->nodes_refine; i++)
            {
                fsX->ncoor_refine[i] = FSLib->rangeX_begin + fsX->dx_refine * i;
            }
            // y-direction
            for(j = 0; j < fsY->nodes_refine; j++)
            {
                fsY->ncoor_refine[j] = FSLib->rangeY_begin + fsY->dx_refine * j;
            }           
        }
    }

    if(1 == FSLib->non_uniform_grid)
    {
        // create a global coordinate for LaMEM
        if(0 == FSLib->fs2D)
        {
            if(2 == mode || 0 == step_fs)
            {
                FSLib->ncoor_ori_x = (PetscScalar *) malloc( (FSLib->nx_fs + 1) * sizeof(PetscScalar));   
                FSLib->ncoor_ori_y = (PetscScalar *) malloc( (FSLib->ny_fs + 1) * sizeof(PetscScalar));  
            }

            // x-direction
            for(i = 0; i < FSLib->msx_fs.nsegs; i++)
            {   
                // generate nodal coordinates for the local part of the segment
                ierr = FastScapeCreateGlobalGrid(FSLib->ncoor_ori_x, FSLib->msx_fs, i, scal); CHKERRQ(ierr);
            } 
            // last nodes
            FSLib->ncoor_ori_x[FSLib->nx_fs] = FSLib->rangeX_end;

            // y-direction
            for(i = 0; i < FSLib->msy_fs.nsegs; i++)
            {   
                // generate nodal coordinates for the local part of the segment
                ierr = FastScapeCreateGlobalGrid(FSLib->ncoor_ori_y, FSLib->msy_fs, i, scal); CHKERRQ(ierr);
            } 
            // last nodes
            FSLib->ncoor_ori_y[FSLib->ny_fs] = FSLib->rangeY_end;   
        }
        else
        {
            if(1 == FSLib->extendedX)
            {
                if(2 == mode || 0 == step_fs)
                {
                    FSLib->ncoor_ori_y = (PetscScalar *) malloc( (FSLib->ny_fs + 1) * sizeof(PetscScalar));
                }

                // y-direction
                for(i = 0; i < FSLib->msy_fs.nsegs; i++)
                {   
                    // generate nodal coordinates for the local part of the segment
                    ierr = FastScapeCreateGlobalGrid(FSLib->ncoor_ori_y, FSLib->msy_fs, i, scal); CHKERRQ(ierr);
                } 
                // last nodes
                FSLib->ncoor_ori_y[FSLib->ny_fs] = FSLib->rangeY_end;                  
            }
            else
            {
                if(2 == mode || 0 == step_fs)
                {
                    FSLib->ncoor_ori_x = (PetscScalar *) malloc( (FSLib->nx_fs + 1) * sizeof(PetscScalar)); 
                }

                // x-direction
                for(i = 0; i < FSLib->msx_fs.nsegs; i++)
                {   
                    // generate nodal coordinates for the local part of the segment
                    ierr = FastScapeCreateGlobalGrid(FSLib->ncoor_ori_x, FSLib->msx_fs, i, scal); CHKERRQ(ierr);
                } 
                // last nodes
                FSLib->ncoor_ori_x[FSLib->nx_fs] = FSLib->rangeX_end;                
            }
        }
    }

    PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode ScalingFastScapeCreate(Scaling *scal)
{
    PetscFunctionBeginUser;
    PetscScalar km, yr, m;

    // read unit values
    km     = 1e3;
    m      = 1e2;
    yr     = 3600.0*24.0*365.25;   

    if( _SI_ == scal->utype)
    {
        // s (LaMEM) -> yr (FastScape)
        scal->time_fs             = 1/yr;                      sprintf(scal->lbl_time_fs,         "[yr]");  
        // m (LaMEM) -> m (FastSCape)
        scal->length_fs           = 1.0;                       sprintf(scal->lbl_length_fs,       "[m]");   
        // m/s (LaMEM) -> m/yr (FastScape)
        scal->velocity_fs         = 1/yr;                      sprintf(scal->lbl_velocity_fs,     "[m/yr]"); 

        // output 
        // m^2 --> m^2
        scal->area_fs             = 1.0;                       sprintf(scal->lbl_area_fs,            "[m^2]"); 
        // m/yr --> m/yr
        scal->rate                = 1.0;                       sprintf(scal->lbl_rate,            "[m/yr]"); 
        scal->degree              = 1.0;                       sprintf(scal->lbl_degree,          "[°]");   
    }
    else if( _GEO_ == scal->utype)
    {
        // Myr (LaMEM) -> yr (FastScape)
        scal->time_fs             = 1e6;                       sprintf(scal->lbl_time_fs,         "[yr]");  
        // km (LaMEM) -> m (FastScape) 
        scal->length_fs           = km;                        sprintf(scal->lbl_length_fs,       "[m]");   
        // cm/yr (LaMEM) -> m/yr (FastScape)
        scal->velocity_fs         = m;                         sprintf(scal->lbl_velocity_fs,     "[m/yr]"); 

        // output 
        // m^2 --> km^2
        scal->area_fs             = km * km;                 sprintf(scal->lbl_area_fs,         "[km^2]"); 
        // m/yr --> km/yr
        scal->rate                = 1.0 * km;                    sprintf(scal->lbl_rate,            "[km/yr]"); 
        scal->degree              = 1.0;                       sprintf(scal->lbl_degree,          "[°]");   
    }
    else
    {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect unit type for FastScape");
    }

    PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PVSurfFastScapeCreate(FastScapeLib *FSLib, FB *fb)
{

    char filename[_str_len_];

    PetscErrorCode ierr;
    PetscFunctionBeginUser;

    // initialize
    // topography & pvd
    FSLib->outsurf_fs            = 1;
    FSLib->outpvd_fs             = 1;
    FSLib->out_topofs            = 1;

    // surface parameter
    FSLib->out_silt_fraction   = 0;
    FSLib->out_basement        = 0;
    FSLib->out_total_erosion   = 0;
    FSLib->out_drainage_area   = 0;
    FSLib->out_erosion_rate    = 0;
    FSLib->out_slope           = 0;
    FSLib->out_curvature       = 0;
    FSLib->out_chi             = 0;
    FSLib->out_catchment       = 0;
    FSLib->out_lake_depth      = 0;

    // check activation
    ierr = getIntParam   (fb, _OPTIONAL_, "out_surf_fs",            &FSLib->outsurf_fs,        1, 1); CHKERRQ(ierr); 
       
    // read
    ierr = getStringParam(fb, _OPTIONAL_, "out_file_name",          filename,               "output"); CHKERRQ(ierr);
    // pvd
    ierr = getIntParam   (fb, _OPTIONAL_, "out_fs_pvd",             &FSLib->outpvd_fs,         1, 1); CHKERRQ(ierr); 

    // topography
    ierr = getIntParam   (fb, _OPTIONAL_, "out_surf_topofs",        &FSLib->out_topofs,            1, 1); CHKERRQ(ierr);

    // surface parameter
    ierr = getIntParam   (fb, _OPTIONAL_, "out_surf_silt_fraction", &FSLib->out_silt_fraction,     1, 1); CHKERRQ(ierr);    
    ierr = getIntParam   (fb, _OPTIONAL_, "out_surf_basement",      &FSLib->out_basement,          1, 1); CHKERRQ(ierr);  
    ierr = getIntParam   (fb, _OPTIONAL_, "out_surf_total_erosion", &FSLib->out_total_erosion,     1, 1); CHKERRQ(ierr);  
    ierr = getIntParam   (fb, _OPTIONAL_, "out_surf_drainage_area", &FSLib->out_drainage_area,     1, 1); CHKERRQ(ierr);  
    ierr = getIntParam   (fb, _OPTIONAL_, "out_surf_erosion_rate",  &FSLib->out_erosion_rate,      1, 1); CHKERRQ(ierr);  
    ierr = getIntParam   (fb, _OPTIONAL_, "out_surf_slope",         &FSLib->out_slope,             1, 1); CHKERRQ(ierr);  
    ierr = getIntParam   (fb, _OPTIONAL_, "out_surf_curvature",     &FSLib->out_curvature,         1, 1); CHKERRQ(ierr);  
    ierr = getIntParam   (fb, _OPTIONAL_, "out_surf_chi",           &FSLib->out_chi,               1, 1); CHKERRQ(ierr);  
    ierr = getIntParam   (fb, _OPTIONAL_, "out_surf_catchment",     &FSLib->out_catchment,         1, 1); CHKERRQ(ierr);  
    ierr = getIntParam   (fb, _OPTIONAL_, "out_surf_lake_depth",    &FSLib->out_lake_depth,        1, 1); CHKERRQ(ierr);  

    // print summary
    PetscPrintf(PETSC_COMM_WORLD, "FastScape output parameters:\n");
    PetscPrintf(PETSC_COMM_WORLD, "   Write .pvd fs file         : %s \n", FSLib->outpvd_fs    ? "yes" : "no");
    if(FSLib->outpvd_fs)
    {
        if(FSLib->out_topofs)        PetscPrintf(PETSC_COMM_WORLD, "     Topo                       @ \n");
        if(FSLib->out_silt_fraction) PetscPrintf(PETSC_COMM_WORLD, "     silt_fraction              @ \n");
        if(FSLib->out_basement)      PetscPrintf(PETSC_COMM_WORLD, "     basement                   @ \n");
        if(FSLib->out_drainage_area) PetscPrintf(PETSC_COMM_WORLD, "     drainage_area              @ \n");
        if(FSLib->out_erosion_rate)  PetscPrintf(PETSC_COMM_WORLD, "     erosion_rate               @ \n");
        if(FSLib->out_slope)         PetscPrintf(PETSC_COMM_WORLD, "     slope                      @ \n");
        if(FSLib->out_curvature)     PetscPrintf(PETSC_COMM_WORLD, "     curvature                  @ \n");
        if(FSLib->out_chi)           PetscPrintf(PETSC_COMM_WORLD, "     chi                        @ \n");
        if(FSLib->out_catchment)     PetscPrintf(PETSC_COMM_WORLD, "     catchment                  @ \n");
        if(FSLib->out_lake_depth)    PetscPrintf(PETSC_COMM_WORLD, "     lake_depth                 @ \n");    
    }
    PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------\n");

    // set file name
    sprintf(FSLib->outfile_fs,        "%s_fs",        filename);

    PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FastScapeCopyVelocity(FastScapeLib *FSLib)
{

    PetscErrorCode ierr;
    JacRes *jr;
    FreeSurf *surf;
    
    surf = FSLib->surf;
    jr = surf->jr;

    ierr = FreeSurfGetVelComp(surf, &InterpXFaceCorner, jr->lvx, surf->vx); CHKERRQ(ierr);
    ierr = FreeSurfGetVelComp(surf, &InterpYFaceCorner, jr->lvy, surf->vy); CHKERRQ(ierr);
    ierr = FreeSurfGetVelComp(surf, &InterpZFaceCorner, jr->lvz, surf->vz); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FastScapeCreateGlobalGrid(PetscScalar *ncoor, MeshSeg1DFS ms_fs, PetscInt iseg, Scaling *scal)
{
    PetscInt    i, M, sum = 0, istart;
    PetscScalar xstart, xclose, bias, avgSz, begSz, endSz, dx;

    PetscFunctionBeginUser;

    // start index
    istart = ms_fs.istart[iseg];

    // total number of cells
    M = ms_fs.istart[iseg+1] - ms_fs.istart[iseg];

    // starting & closing coordinates
    xstart = ms_fs.xstart[iseg] * scal->length;
    xclose = ms_fs.xstart[iseg+1] * scal->length;

    // bias (last to first cell size ratio > 1 -> growing)
    bias = ms_fs.biases[iseg];

    // average cell size
    avgSz = (xclose - xstart)/(PetscScalar)M;

    // uniform case
    if(bias == 1.0)
    {
        // generate coordinates of local nodes
        for(i = istart; i < M + istart + 1; i++)
        {   
            ncoor[i] = xstart + (PetscScalar)( (i - istart) * avgSz );
        }
    }
    // non-uniform case
    else
    {
        // cell size limits
        begSz = 2.0*avgSz/(1.0 + bias);
        endSz = bias*begSz;

        // cell size increment (negative for bias < 1)
        dx = (endSz - begSz)/(PetscScalar)(M-1);
		
		// generate coordinates of local nodes
		for(i = istart; i < istart + M; i++)
		{
			ncoor[i] = xstart + (i - istart) * begSz + dx * sum ;
            sum += i - istart;
		}
    }

    // override last node coordinate
    ncoor[ ms_fs.istart[iseg+1] ] = xclose;
    
    PetscFunctionReturn(0); 
}
//---------------------------------------------------------------------------
PetscErrorCode InterpolationFor3DNonUniformGrid(FastScapeLib *FSLib, PetscScalar *value, PetscInt mode)
{
    FSGrid    *fsX;
    FSGrid    *fsY;
    fsX     = &FSLib->fsX;
    fsY     = &FSLib->fsY;

    PetscErrorCode ierr;       
    PetscInt i, j, m, n, mm, nn, ind, ind_a, ind_b, ind_c, ind_d, p, q, find_m = 0, find_n = 0, find_indx = 0, find_indy = 0;
    PetscScalar dx, dy, wtx, wty, x1, x2, y1, y2, x_coor, y_coor;
    PetscScalar *value_save;

    if(1 == mode)
    {
        ierr = VecGetArray(FSLib->gtopo_nug, &value_save);       CHKERRQ(ierr);    
    }
    else if(2 == mode)
    {
        ierr = VecGetArray(FSLib->vx_nug, &value_save);          CHKERRQ(ierr); 
    }
    else if(3 == mode)
    {
        ierr = VecGetArray(FSLib->vy_nug, &value_save);          CHKERRQ(ierr); 
    }   
    else if(4 == mode)
    { 
        ierr = VecGetArray(FSLib->vz_nug, &value_save);          CHKERRQ(ierr); 
    }

    for(j = 0; j < fsY->nodes; j++)
    {
        for(i = 0; i < fsX->nodes; i++)
        {     
            x_coor = fsX->ncoor[i];
            y_coor = fsY->ncoor[j];   

            // get nearest four index
            for(p = 0; p < FSLib->nx_fs; )
            {
                x1 = FSLib->ncoor_ori_x[p];
                x2 = FSLib->ncoor_ori_x[p + 1];

                if( (x1 <= x_coor) && (x2 >= x_coor) )
                {
                    m      = p;
                    mm     = p + 1;
                    if( FSLib->nx_fs == mm ) 
                    {
                        mm -= 1;
                    }
                        
                    find_m = 1;                   
                    break;
                }   
                 
                find_indx = floor( (x_coor - x1) / FSLib->msx_fs.max_spacing );
                if( 0 < find_indx)
                {
                    p += find_indx;
                }
                else
                {
                    p++;
                }
            }

            if( 1 == find_m )
            {
                for(q = 0; q < FSLib->ny_fs; )
                {
                    y1 = FSLib->ncoor_ori_y[q];
                    y2 = FSLib->ncoor_ori_y[q + 1];
                        
                    if( (y1 <= y_coor) && (y2 >= y_coor) )
                    {
                        n      = q;
                        nn     = q + 1;
                        if( FSLib->ny_fs == nn ) 
                        {
                            nn -= 1;
                        }
                        find_n = 1;                      
                        break;
                    }               
  
                    find_indy = floor( (y_coor - y1) / FSLib->msy_fs.max_spacing );
                    if( 0 < find_indy)
                    {
                        q += find_indy;
                    }
                    else
                    {
                        q++;
                    }
                }
                
                if( 1 == find_n )
                {
                    // interpolate
                    ind = j * fsX->nodes + i;    
                    ind_a = n * FSLib->nx_fs + m;
                    ind_b = n * FSLib->nx_fs + mm;
                    ind_c = nn * FSLib->nx_fs + m;
                    ind_d = nn * FSLib->nx_fs + mm;                             
                    // bilinear interpolation
                    dx    = x2 - x1;
                    dy    = y2 - y1;
                    wtx   = x_coor - x1;
                    wty   = y_coor - y1;
                     
                    // boundary (bottom or right)
                    if(m == mm)
                    {
                        if(n != nn)
                        {
                            value_save[ind] = ( (dy - wty) / dy )  * value[ind_a]+ \
                                                ( wty / dy )  * value[ind_c];  
                        }                     
                    }
                    else if(n == nn)
                    {
                        if(m != mm)
                        {
                            value_save[ind] = ( (dx - wtx) / dx ) * value[ind_a]+ \
                                                ( wtx / dx ) * value[ind_b];
                        }
                        else
                        {
                            value_save[ind] = value[ind_a];
                        }
                    }
                    else
                    {
                        value_save[ind] = ( (dy - wty) / dy ) * ( (dx - wtx) / dx ) * value[ind_a]+ \
                                            ( (dy - wty) / dy ) * ( wtx / dx ) * value[ind_b] + \
                                            ( wty / dy ) * ( (dx - wtx) / dx ) * value[ind_c] + \
                                            ( wty / dy ) * ( wtx / dx ) * value[ind_d]; 
                    }
                } 
            }
                        
            find_m = 0;
            find_n = 0;
                
        }
    }
    
    if(1 == mode)
    {
        ierr = VecRestoreArray(FSLib->gtopo_nug, &value_save);       CHKERRQ(ierr);    
    }
    else if(2 == mode)
    {
        ierr = VecRestoreArray(FSLib->vx_nug, &value_save);          CHKERRQ(ierr); 
    }
    else if(3 == mode)
    {
        ierr = VecRestoreArray(FSLib->vy_nug, &value_save);          CHKERRQ(ierr); 
    }   
    else if(4 == mode)
    { 
        ierr = VecRestoreArray(FSLib->vz_nug, &value_save);          CHKERRQ(ierr); 
    }

    PetscFunctionReturn(0);    
}
//---------------------------------------------------------------------------
PetscErrorCode InterpolationFor2DNonUniformGrid(FastScapeLib *FSLib, PetscScalar *value_in, PetscScalar *value_out)
{
    FSGrid    *fsX;
    FSGrid    *fsY;
    fsX     = &FSLib->fsX;
    fsY     = &FSLib->fsY;
    
    PetscInt i, j, m, n, mm, nn, ind, ind_a, ind_b, p, q, find_m = 0, find_n = 0, find_indx = 0, find_indy = 0;
    PetscScalar dx, dy, wtx, wty, x1, x2, y1, y2, x_coor, y_coor;

    // interpolate after calculating mean value
    if( 1 == FSLib->extendedX)
    {
        // get nearest 2 index
       for(j = 0; j < fsY->nodes_extend; j++)
       {
            y_coor = fsY->ncoor_extend[j];  

            for(q = 0; q < FSLib->ny_fs; )
            {
                y1 = FSLib->ncoor_ori_y[q];
                y2 = FSLib->ncoor_ori_y[q + 1];
                
                if( (y1 <= y_coor) && (y2 >= y_coor) )
                {
                    n      = q;
                    nn     = q + 1;
                    if( FSLib->ny_fs == nn ) 
                    {
                        nn -= 1;
                    }
                    find_n = 1;                      
                    break;
                }               

                find_indy = floor( (y_coor - y1) / FSLib->msy_fs.max_spacing );
                if( 0 < find_indy)
                {
                    q += find_indy;
                }
                else
                {
                    q++;
                }
            }
            
            if( 1 == find_n )
            {
                // interpolate
                ind = j;    
                ind_a = n;
                ind_b = nn;  

                // linear interpolation
                dy    = y2 - y1;
                wty   = y_coor - y1;

                if(n == nn)
                {
                    value_out[ind] = value_in[ind_a];
                }
                else
                {
                    value_out[ind] = ( (dy - wty) / dy ) * value_in[ind_a]+ ( wty / dy ) * value_in[ind_b];
                }

            }                 
            find_n = 0;
        }
    }
    else
    {
        // get nearest 2 index
        for(i = 0; i < fsX->nodes_extend; i++)
        {
            x_coor = fsX->ncoor_extend[i];

            // get nearest 2 index
            for(p = 0; p < FSLib->nx_fs; )
            {
                x1 = FSLib->ncoor_ori_x[p];
                x2 = FSLib->ncoor_ori_x[p + 1];

                if( (x1 <= x_coor) && (x2 >= x_coor) )
                {
                    m      = p;
                    mm     = p + 1;
                    if( FSLib->nx_fs == mm ) 
                    {
                        mm -= 1;
                    }
                        
                    find_m = 1;                   
                    break;
                }   
                   
                find_indx = floor( (x_coor - x1) / FSLib->msx_fs.max_spacing );
                if( 0 < find_indx)
                {
                    p += find_indx;
                }
                else
                {
                    p++;
                }
            }

            if( 1 == find_m )
            {
                // interpolate
                ind = i;    
                ind_a = m;
                ind_b = mm;                          
                // bilinear interpolation
                dx    = x2 - x1;
                wtx   = x_coor - x1;
                
                if(m == mm)
                {
                    value_out[ind] = value_in[ind_a];
                }
                else
                {
                    value_out[ind] = ( (dx - wtx) / dx ) * value_in[ind_a]+ ( wtx / dx ) * value_in[ind_b];
                }
            } 

            find_m = 0;
        }
    }
  
    PetscFunctionReturn(0);    
}
//---------------------------------------------------------------------------
PetscErrorCode GatherVariableFromLaMEM(FastScapeLib *FSLib, PetscScalar *topo_alloc, PetscScalar *vx_alloc, PetscScalar *vy_alloc, PetscScalar *vz_alloc, PetscInt step_fs)
{
    PetscErrorCode ierr;
    PetscInt tnodes, L, sx, sy, sz, nx, ny, nz, tproc, cproc, rankZ_id, i, j, k, ind;
    PetscScalar vz_f, vz_f2;
    PetscScalar ***topo;
    PetscScalar ***vz, ***vx, ***vy;
    PetscScalar ***vz_collect, ***vx_collect, ***vy_collect;
    PetscScalar *topo_collect = PETSC_NULL;
    PetscScalar *vz_fs = PETSC_NULL, *vx_fs = PETSC_NULL, *vy_fs = PETSC_NULL;
  
    FDSTAG      *fs;
    FreeSurf    *surf;

    surf   = FSLib->surf;
    fs     = surf->jr->fs;

    tnodes = FSLib->nx_fs * FSLib->ny_fs;

    // Gather topography and velocity
    // topography
    ierr = DMDAVecGetArray(surf->DA_SURF, surf->gtopo,  &topo);  CHKERRQ(ierr);
    
    if(0 == step_fs)
    {
        ierr = VecScatterCreateToZero(surf->gtopo, &FSLib->ctx, &FSLib->gtopo_collect);                        CHKERRQ(ierr); 
        ierr = VecScatterBegin(FSLib->ctx, surf->gtopo, FSLib->gtopo_collect, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = VecScatterEnd(FSLib->ctx, surf->gtopo, FSLib->gtopo_collect, INSERT_VALUES, SCATTER_FORWARD);   CHKERRQ(ierr);
    }

    // velocity
    // vx & vy
    ierr = DMDAVecGetArray(surf->DA_SURF, surf->vx, &vx);                      CHKERRQ(ierr);
    ierr = DMDAVecGetArray(surf->DA_SURF, FSLib->vx_collect, &vx_collect);     CHKERRQ(ierr);
    ierr = DMDAVecGetArray(surf->DA_SURF, surf->vy, &vy);                      CHKERRQ(ierr);
    ierr = DMDAVecGetArray(surf->DA_SURF, FSLib->vy_collect, &vy_collect);     CHKERRQ(ierr);
    // vz
    L    = (PetscInt)fs->dsz.rank;
    ierr = DMDAGetCorners(fs->DA_COR, &sx, &sy, &sz, &nx, &ny, &nz);           CHKERRQ(ierr);
    ierr = DMDAVecGetArray(surf->DA_SURF, surf->vz, &vz);                      CHKERRQ(ierr);
    ierr = DMDAVecGetArray(surf->DA_SURF, FSLib->vz_collect, &vz_collect);     CHKERRQ(ierr);

    // load local value to global varibles
    START_PLANE_LOOP
    {
        vx_collect[L][j][i] = vx[L][j][i]; // (km)
        vy_collect[L][j][i] = vy[L][j][i]; // (km)
        vz_collect[L][j][i] = vz[L][j][i]; // (km)        
    }
    END_PLANE_LOOP  

    // note: if the vec is a local vec, it will create a local output, which isn't correct
    // A global vector is needed
    ierr = VecScatterCreateToZero(FSLib->vx_collect, &FSLib->ctx, &FSLib->vx_fs);                        CHKERRQ(ierr); 
    ierr = VecScatterBegin(FSLib->ctx, FSLib->vx_collect, FSLib->vx_fs, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(FSLib->ctx, FSLib->vx_collect, FSLib->vx_fs, INSERT_VALUES, SCATTER_FORWARD);   CHKERRQ(ierr);

    ierr = VecScatterCreateToZero(FSLib->vy_collect, &FSLib->ctx, &FSLib->vy_fs);                        CHKERRQ(ierr); 
    ierr = VecScatterBegin(FSLib->ctx, FSLib->vy_collect, FSLib->vy_fs, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(FSLib->ctx, FSLib->vy_collect, FSLib->vy_fs, INSERT_VALUES, SCATTER_FORWARD);   CHKERRQ(ierr);
    
    ierr = VecScatterCreateToZero(FSLib->vz_collect, &FSLib->ctx, &FSLib->vz_fs);                            CHKERRQ(ierr); 
    ierr = VecScatterBegin(FSLib->ctx, FSLib->vz_collect, FSLib->vz_fs, INSERT_VALUES, SCATTER_FORWARD);     CHKERRQ(ierr);
    ierr = VecScatterEnd(FSLib->ctx, FSLib->vz_collect, FSLib->vz_fs, INSERT_VALUES, SCATTER_FORWARD);       CHKERRQ(ierr);

    // reallocate
    tproc = fs->dsx.nproc * fs->dsy.nproc * fs->dsz.nproc;

    PetscInt para_info[7] = {sx, sy, sz, nx, ny, nz, 0};

    MPI_Comm_rank(PETSC_COMM_WORLD, &cproc);

    vz_f = vz[L][sy][sx];
    vz_f2 = vz[L][sy][sx+1];

    // send message to processor 0
    if(!ISRankZero(PETSC_COMM_WORLD))
    {
        if(0 == fs->dsz.rank)
        {
            // current order
            if(1 == fs->dsx.nproc)
            {
                rankZ_id = fs->dsy.rank;
            }
            else
            { 
                rankZ_id = fs->dsy.nproc * fs->dsy.rank + fs->dsx.rank;
            }
        }
        else
        {
            rankZ_id = NO_NEED;
        }

        para_info[6] = rankZ_id;

        ierr = MPI_Send(para_info, 7, MPIU_INT, 0, 0, PETSC_COMM_WORLD); CHKERRQ(ierr); 
        ierr = MPI_Send(&vz_f, 1, MPIU_SCALAR, 0, 1, PETSC_COMM_WORLD);  CHKERRQ(ierr);
        ierr = MPI_Send(&vz_f2, 1, MPIU_SCALAR, 0, 2, PETSC_COMM_WORLD); CHKERRQ(ierr);
    }

    // recieve message in rank 0
    if(ISRankZero(PETSC_COMM_WORLD))
    {
        PetscInt    countI = 0, ind_alloc = 0, countJ;
        PetscInt    **para_info_rec = PETSC_NULL;
        PetscScalar *vz_first = PETSC_NULL, *vz_first2 = PETSC_NULL;

        vz_first     = (PetscScalar *)malloc(tproc * sizeof(PetscScalar)); 
        vz_first2    = (PetscScalar *)malloc(tproc * sizeof(PetscScalar));

        para_info_rec = (PetscInt **)malloc(tproc * sizeof(PetscInt*));
        for(i = 0; i < tproc; i++) para_info_rec[i] = (PetscInt *) malloc(7 * sizeof(PetscInt));

        for(i = 0; i < tproc; i++)
        {
            if(0 == i)
            {
                for(j = 0; j < 7; j++)
                {
                    para_info_rec[i][j] = para_info[j];
                }

                vz_first[i] = vz[L][0][0];
                vz_first2[i] = vz[L][0][1];
            }
            else
            {
                ierr = MPI_Recv(para_info_rec[i], 7, MPIU_INT, i, 0, PETSC_COMM_WORLD, MPI_STATUS_IGNORE); CHKERRQ(ierr);
                ierr = MPI_Recv(&vz_first[i], 1, MPIU_SCALAR, i, 1, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);  CHKERRQ(ierr);
                ierr = MPI_Recv(&vz_first2[i], 1, MPIU_SCALAR, i, 2, PETSC_COMM_WORLD, MPI_STATUS_IGNORE); CHKERRQ(ierr);
            }
        }

        // reallocate
        // All the data of a process is put into a new array in a whole, and the index may not be correct
        if( 0 == step_fs )
        {
            ierr = VecGetArray(FSLib->gtopo_collect, &topo_collect);  CHKERRQ(ierr);
        }

        ierr = VecGetArray(FSLib->vx_fs,  &vx_fs);  CHKERRQ(ierr);
        ierr = VecGetArray(FSLib->vy_fs,  &vy_fs);  CHKERRQ(ierr);        
        ierr = VecGetArray(FSLib->vz_fs, &vz_fs);       CHKERRQ(ierr);

        do
        {
            for(i = 0; i < tproc; i++)
            {
                if( (vz_fs[ind_alloc] == vz_first[i]) && (vz_fs[ind_alloc+1] == vz_first2[i]) && (para_info_rec[i][6] != NO_NEED) )
                {

                    countJ = 0;          

                    for(j = para_info_rec[i][1]; j < para_info_rec[i][1] + para_info_rec[i][4]; j++)
                    {
                        for(k = para_info_rec[i][0]; k < para_info_rec[i][0] + para_info_rec[i][3]; k++)
                        {
                            ind = j * FSLib->nx_fs + k;

                            // topography
                            if( 0 == step_fs ) topo_alloc[ind] = topo_collect[countI]; 
                            
                            // vx, vy, vz
                            vx_alloc[ind] = vx_fs[countI];    
                            vy_alloc[ind] = vy_fs[countI];                                                         
                            vz_alloc[ind] = vz_fs[countI];

                            countI++;
                            countJ++;
                        }
                    }

                    ind_alloc += countJ; 
                }
            }
        } while (ind_alloc < tnodes);
        

        free(vz_first);
        free(vz_first2);
        free(para_info_rec);
 
        vz_first = PETSC_NULL;
        vz_first2 = PETSC_NULL;
        para_info_rec = PETSC_NULL;

        if( 0 == step_fs )
        {
            ierr = VecRestoreArray(FSLib->gtopo_collect, &topo_collect);      CHKERRQ(ierr);
        }

        ierr = VecRestoreArray(FSLib->vz_fs, &vz_fs);       CHKERRQ(ierr);
        ierr = VecRestoreArray(FSLib->vx_fs, &vx_fs);       CHKERRQ(ierr);
        ierr = VecRestoreArray(FSLib->vy_fs, &vy_fs);       CHKERRQ(ierr);
        

    }

    ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->gtopo, &topo);              CHKERRQ(ierr);

    ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->vx, &vx);                     CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(surf->DA_SURF, FSLib->vx_collect, &vx_collect);    CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->vy, &vy);                     CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(surf->DA_SURF, FSLib->vy_collect, &vy_collect);    CHKERRQ(ierr);         
    ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->vz, &vz);                   CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(surf->DA_SURF, FSLib->vz_collect, &vz_collect);  CHKERRQ(ierr);    
    
    ierr = VecScatterDestroy(&FSLib->ctx); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FastScapeStretchGrid(FastScapeLib *FSLib)
{
    PetscErrorCode ierr;
    PetscScalar Exx, Eyy, Ezz, Rxx, Ryy, Rzz, step;
    PetscInt i, j;
    Scaling  *scal;
    JacRes   *jr;

    jr    = FSLib->jr;
    scal  = FSLib->scal;
    step  = jr->ts->dt;

    // update model range
    // range X, Y, Z   
    ierr  = BCGetBGStrainRates(jr->bc, &Exx, &Eyy, &Ezz, NULL, NULL, NULL, &Rxx, &Ryy, &Rzz); CHKERRQ(ierr);

    FSLib->rangeX_begin += (step * Exx * (FSLib->rangeX_begin / scal->length - Rxx) ) * scal->length; // km (LaMEM) in GEO
    FSLib->rangeX_end   += (step * Exx * (FSLib->rangeX_end / scal->length - Rxx) ) * scal->length;
    FSLib->rangeY_begin += (step * Eyy * (FSLib->rangeY_begin / scal->length - Ryy) ) * scal->length;
    FSLib->rangeY_end   += (step * Eyy * (FSLib->rangeY_end / scal->length - Ryy) ) * scal->length;
    FSLib->rangeZ_begin += (step * Ezz * (FSLib->rangeZ_begin / scal->length - Rzz) ) * scal->length;   
    FSLib->rangeZ_end   += (step * Ezz * (FSLib->rangeZ_end / scal->length - Rzz) ) * scal->length;

    FSLib->rangeX       = (FSLib->rangeX_end - FSLib->rangeX_begin) * scal->length_fs; //(km) in LaMEM to (m) in FastScape (GEO)
    FSLib->rangeY       = (FSLib->rangeY_end - FSLib->rangeY_begin) * scal->length_fs;

    if(1 == FSLib->fs2D)
    {
        // 2D No Refine
        // extend grid in FastScape
        // extend in rangeX
        if(1 == FSLib->extendedX)
        {
            if(Eyy)
            {
                FSLib->extendedXRange += (step * Exx * (FSLib->extendedXRange / scal->length - Rxx)) * scal->length * scal->length_fs;
                FSLib->extendedYRange = FSLib->rangeY;
            }
            else
            {
                FSLib->extendedXRange += FSLib->rangeX_begin * scal->length_fs;
                FSLib->extendedYRange = FSLib->rangeY;                
            }
        }
        else
        {
            if(Exx)
            {
                FSLib->extendedXRange = FSLib->rangeX;
                FSLib->extendedYRange += (step * Eyy * (FSLib->extendedYRange / scal->length - Ryy)) * scal->length * scal->length_fs;
            }
            else
            {
                FSLib->extendedXRange = FSLib->rangeX;
                FSLib->extendedYRange += FSLib->rangeY_begin * scal->length_fs;                
            }
        }
    }

    if(1 == FSLib->non_uniform_grid)
    { 
        if(0 == FSLib->fs2D)
        {
            // x-direction
            for(i = 0; i < FSLib->msx_fs.nsegs + 1; i++)
            {
                FSLib->msx_fs.xstart[i] += step * Exx * (FSLib->msx_fs.xstart[i] - Rxx);
            }

            // y-direction
            for(j = 0; j < FSLib->msy_fs.nsegs + 1; j++)
            {
                FSLib->msy_fs.xstart[j] += step * Eyy * (FSLib->msy_fs.xstart[j] - Ryy);
            }
        }   
        else
        {
            if(1 == FSLib->extendedX)
            {
                // y-direction
                for(j = 0; j < FSLib->msy_fs.nsegs + 1; j++)
                {
                    FSLib->msy_fs.xstart[j] += step * Eyy * (FSLib->msy_fs.xstart[j] - Ryy);
                }
            }
            else
            {
                // x-direction
                for(i = 0; i < FSLib->msx_fs.nsegs + 1; i++)
                {
                    FSLib->msx_fs.xstart[i] += step * Exx * (FSLib->msx_fs.xstart[i] - Rxx);
                }
            }
        }
    }

    // update coordinate
    ierr = FastScapeCreateSurfaceGrid(FSLib, 0); CHKERRQ(ierr);

    PetscFunctionReturn(0); 
}
//---------------------------------------------------------------------------
PetscErrorCode BilinearInterpolate(FastScapeLib *FSLib, PetscScalar *data, PetscScalar *data_refine, 
    PetscScalar *data_refine_pass, Scaling *scal, PetscInt corMode, PetscInt nx_refine, PetscInt ny_refine)
{
    //PetscInt interpolationMode = 1;
    // linear interpolation
    // refine = 1; no nodes; 2, add a node between two original nearest nodes; 3, add two nodes between two nearest original nodes;
    // corMode: 1 -- topography, (km) in LaMEM to (m) in FastScape; 2 -- velocity, (cm/yr) in LaMEM to (m/yr) in FastScape
    PetscInt i, j, ind, ind_a, ind_b, ind_aa, ind_bb, countX1, countX2, countY, tnodes_refine;
    PetscScalar distance_ax, distance_ay;

    tnodes_refine = nx_refine * ny_refine;
    countX1       = 0;
    countX2       = 0;
    countY        = 0;
    ind_a         = countX1 * FSLib->refine;
    ind_b         = (countX1 + 1) * FSLib->refine;
    ind_aa        = countX2;
    ind_bb        = countX2 + 1;

    // interpolate in x-direction
    for(j = 0; j < ny_refine; j += FSLib->refine) 
    {   
        for(i = 0; i < nx_refine; i++) 
        {
            ind = j * nx_refine + i;

            if(0 == ind % nx_refine) 
            { 
                countX1 = 0;   

                if(0 == ind)
                {
                    if(1 == corMode)
                    {

                            data_refine[ind_a] = data[ind_aa] * scal->length * scal->length_fs;
                            data_refine[ind_b] = data[ind_bb] * scal->length * scal->length_fs;
                                            
                    }
                    else if(2 == corMode)
                    {
                        data_refine[ind_a] = data[ind_aa] * scal->velocity / scal->velocity_fs;
                        data_refine[ind_b] = data[ind_bb] * scal->velocity / scal->velocity_fs; 
                    }
                    else
                    {
                        data_refine[ind_a] = data[ind_aa];
                        data_refine[ind_b] = data[ind_bb];                          
                    }                
                }
                else
                {
                    ind_a  = countX1 * FSLib->refine + j * nx_refine;
                    ind_b  = (countX1+1) * FSLib->refine + j * nx_refine; 
                    
                    countX2++; 
                    
                    ind_aa = countX2;
                    ind_bb = countX2 + 1;
 
                    if(1 == corMode)
                    {

                            data_refine[ind_a] = data[ind_aa] * scal->length * scal->length_fs;
                            data_refine[ind_b] = data[ind_bb] * scal->length * scal->length_fs;  
                                               
                    }
                    else if(2 == corMode)
                    {
                        data_refine[ind_a] = data[ind_aa] * scal->velocity / scal->velocity_fs;
                        data_refine[ind_b] = data[ind_bb] * scal->velocity / scal->velocity_fs; 
                    }
                    else            
                    {
                        data_refine[ind_a] = data[ind_aa];
                        data_refine[ind_b] = data[ind_bb];                        
                    }        

                }
                
                data_refine_pass[ind_a] = data_refine[ind_a];
                data_refine_pass[ind_b] = data_refine[ind_b]; 
            }
            else
            {
                if(ind == ind_b) 
                {
                    if(ind == tnodes_refine-1) 
                    {
                        if(1 == corMode)
                        {

                                data_refine[ind_b] = data[ind_bb] * scal->length * scal->length_fs;                                 
                            
                                                   
                        }
                        else if(2 == corMode)
                        {
                            data_refine[ind_b] = data[ind_bb] * scal->velocity / scal->velocity_fs; 
                        }
                        else 
                        {
                            data_refine[ind_b] = data[ind_bb];
                        }
                        
                        goto skip;
                    }
                    else
                    {
                        countX1++;
                        countX2++;

                        ind_a  = countX1 * FSLib->refine + j * nx_refine;
                        ind_b  = (countX1+1) * FSLib->refine + j * nx_refine;
                        ind_aa = countX2;
                        ind_bb = countX2 + 1;

                        if(1 == corMode)
                        {

                                data_refine[ind_b] = data[ind_bb] * scal->length * scal->length_fs;  
                                                    
                        }
                        else if(2 == corMode)
                        {
                            data_refine[ind_b] = data[ind_bb] * scal->velocity / scal->velocity_fs; 
                        }
                        else  
                        {
                            data_refine[ind_b] = data[ind_bb];
                        }
                    }

                    data_refine_pass[ind_b]   = data_refine[ind_b];
                }
                else
                {
                    distance_ax           = ((ind - ind_a) % nx_refine) * 1.0 / FSLib->refine;   
                    data_refine[ind]      = data_refine[ind_a] * (1 - distance_ax)  + data_refine[ind_b] * distance_ax;               
                    data_refine_pass[ind] = data_refine[ind];
                }
            }                   
        }
    }
    skip: // printf("interpolation in x-direction done\n");     

    // interpolate in y-direction
    for(j = 0; j < ny_refine; j++) 
    {   
        if(0 == j % FSLib->refine) 
        {
            if( 0 != j) countY += FSLib->refine;    
            continue;
        }

        for(i = 0; i < nx_refine; i++) 
        {

            ind   = j * nx_refine + i;
            ind_a = countY * nx_refine + i;
            ind_b = (countY + FSLib->refine) * nx_refine + i;

            distance_ay           = (j - countY) * 1.0 / FSLib->refine;  
            data_refine[ind]      = data_refine[ind_a] * (1 - distance_ay)  + data_refine[ind_b] * distance_ay;                   
            data_refine_pass[ind] = data_refine[ind]; 
        }         
    }   

    PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode Extended2D(FastScapeLib *FSLib, PetscScalar *data, PetscScalar *data_extended, PetscScalar *data_extend_pass, Scaling *scal, PetscInt corMode)
{
    PetscScalar sum;
    PetscInt i, j, ind; 
    PetscErrorCode ierr;
    PetscScalar *data_aver_ori = PETSC_NULL;
    PetscScalar *data_aver     = PETSC_NULL;

    // extend in x-direction
    if(1 == FSLib->extendedX)
    {
        data_aver_ori = (PetscScalar *)malloc(FSLib->ny_fs           * sizeof(PetscScalar)); 
        data_aver     = (PetscScalar *)malloc(FSLib->extendedYNodes  * sizeof(PetscScalar)); 

        // get mean values
        for(j = 0; j < FSLib->ny_fs; j++)
        {
            sum = 0;

            for(i = 0; i < FSLib->nx_fs; i++)
            {
                ind = j * FSLib->nx_fs + i;

                if(1 == corMode)
                {
                    sum = sum + data[ind] * scal->length * scal->length_fs; // km -> m
                }
                if(2 == corMode)
                {
                    sum = sum + data[ind] * scal->velocity / scal->velocity_fs; // cm/yr -> m/yr
                }
            }

            data_aver_ori[j] = sum / FSLib->nx_fs;
        }        

        if(0 == FSLib->non_uniform_grid)
        {
            for(j = 0; j < FSLib->extendedYNodes; j++)
            {
                for(i = 0; i < FSLib->extendedXNodes; i++)
                {
                    ind = j * FSLib->extendedXNodes + i;

                    data_extended[ind]    = data_aver_ori[j];
                    data_extend_pass[ind] = data_extended[ind];
                }
            }            
        }
        else
        {
            ierr = InterpolationFor2DNonUniformGrid(FSLib, data_aver_ori, data_aver); CHKERRQ(ierr);

            for(j = 0; j < FSLib->extendedYNodes; j++)
            {
                for(i = 0; i < FSLib->extendedXNodes; i++)
                {
                    ind = j * FSLib->extendedXNodes + i;

                    data_extended[ind]    = data_aver[j];
                    data_extend_pass[ind] = data_extended[ind];
                }
            }            
        }
    }
    // extend in y-direction
    else
    {
        data_aver_ori = (PetscScalar *)malloc(FSLib->nx_fs          * sizeof(PetscScalar)); 
        data_aver     = (PetscScalar *)malloc(FSLib->extendedXNodes * sizeof(PetscScalar));  

        for(i = 0; i < FSLib->nx_fs; i++)
        {
            sum = 0;

            for(j = 0; j < FSLib->ny_fs; j++)
            {
                ind = j * FSLib->nx_fs + i;

                if(1 == corMode)
                {
                    sum = sum + data[ind] * scal->length * scal->length_fs; // km -> m
                }
                if(2 == corMode)
                {
                    sum = sum + data[ind] * scal->velocity / scal->velocity_fs; // cm/yr -> m/yr
                }
            }

            data_aver_ori[i] = sum / FSLib->ny_fs;
        }  

        if(0 == FSLib->non_uniform_grid)
        {
            for(i = 0; i < FSLib->extendedXNodes; i++)
            {
                for(j = 0; j < FSLib->extendedYNodes; j++)
                {
                    ind = j * FSLib->extendedXNodes + i;

                    data_extended[ind]    = data_aver_ori[i];
                    data_extend_pass[ind] = data_extended[ind];
                }
            } 
        }
        else
        {
            ierr = InterpolationFor2DNonUniformGrid(FSLib, data_aver_ori, data_aver); CHKERRQ(ierr);
        
            for(i = 0; i < FSLib->extendedXNodes; i++)
            {
                for(j = 0; j < FSLib->extendedYNodes; j++)
                {
                    ind = j * FSLib->extendedXNodes + i;

                    data_extended[ind]    = data_aver[i];
                    data_extend_pass[ind] = data_extended[ind];
                }
            }
        }
    }

    free(data_aver_ori);
    free(data_aver);
    data_aver_ori = PETSC_NULL;
    data_aver = PETSC_NULL;

    PetscFunctionReturn(0);    
}
//---------------------------------------------------------------------------
PetscErrorCode Extended2DRefine(FastScapeLib *FSLib, PetscScalar *data, PetscScalar *data_extended)
{
    // no unit transform
    PetscScalar sum;
    PetscInt i, j, ind; 
    PetscErrorCode ierr;
    PetscScalar *data_aver_ori = PETSC_NULL;
    PetscScalar *data_aver     = PETSC_NULL;

    // extend in x-direction
    if(1 == FSLib->extendedX)
    {
        data_aver_ori = (PetscScalar *)malloc(FSLib->ny_fs           * sizeof(PetscScalar)); 
        data_aver     = (PetscScalar *)malloc(FSLib->extendedYNodes  * sizeof(PetscScalar)); 

        // get mean values
        for(j = 0; j < FSLib->ny_fs; j++)
        {
            sum = 0;

            for(i = 0; i < FSLib->nx_fs; i++)
            {
                ind = j * FSLib->nx_fs + i;

                sum = sum + data[ind]; 
            }

            data_aver_ori[j] = sum / FSLib->nx_fs;
        }        

        if(0 == FSLib->non_uniform_grid)
        {
            for(j = 0; j < FSLib->extendedYNodes; j++)
            {
                for(i = 0; i < FSLib->extendedXNodes; i++)
                {
                    ind = j * FSLib->extendedXNodes + i;

                    data_extended[ind]    = data_aver_ori[j];
                }
            }            
        }
        else
        {
            ierr = InterpolationFor2DNonUniformGrid(FSLib, data_aver_ori, data_aver); CHKERRQ(ierr);

            for(j = 0; j < FSLib->extendedYNodes; j++)
            {
                for(i = 0; i < FSLib->extendedXNodes; i++)
                {
                    ind = j * FSLib->extendedXNodes + i;

                    data_extended[ind]    = data_aver[j];
                }
            }            
        }
    }
    // extend in y-direction
    else
    {
        data_aver_ori = (PetscScalar *)malloc(FSLib->nx_fs          * sizeof(PetscScalar)); 
        data_aver     = (PetscScalar *)malloc(FSLib->extendedXNodes * sizeof(PetscScalar));  

        for(i = 0; i < FSLib->nx_fs; i++)
        {
            sum = 0;

            for(j = 0; j < FSLib->ny_fs; j++)
            {
                ind = j * FSLib->nx_fs + i;

                sum = sum + data[ind]; 
            }

            data_aver_ori[i] = sum / FSLib->ny_fs;
        }  

        if(0 == FSLib->non_uniform_grid)
        {
            for(i = 0; i < FSLib->extendedXNodes; i++)
            {
                for(j = 0; j < FSLib->extendedYNodes; j++)
                {
                    ind = j * FSLib->extendedXNodes + i;

                    data_extended[ind]    = data_aver_ori[i];
                }
            } 
        }
        else
        {
            ierr = InterpolationFor2DNonUniformGrid(FSLib, data_aver_ori, data_aver); CHKERRQ(ierr);
        
            for(i = 0; i < FSLib->extendedXNodes; i++)
            {
                for(j = 0; j < FSLib->extendedYNodes; j++)
                {
                    ind = j * FSLib->extendedXNodes + i;

                    data_extended[ind]    = data_aver[i];
                }
            }
        }
    }

    free(data_aver_ori);
    free(data_aver);
    data_aver_ori = PETSC_NULL;
    data_aver = PETSC_NULL;

    PetscFunctionReturn(0);    
}
//---------------------------------------------------------------------------
PetscErrorCode FastScapeRun(FastScapeLib *FSLib)
{
    /* Unit
    FastScape:
    model range (rangeX, rangeY): m
    timestep: yr
    velocity: m/yr
    topography: m

    LaMEM: when using Geo unit
    model range : km
    timestep after scalling: Myr
    velocity: after scalling: cm/yr
    topography: km
    */

    // Apply surface process to the internal free surface of the model
    // free surface cases only
    FreeSurf *surf;
    surf = FSLib->surf;
    if(!surf->UseFreeSurf) PetscFunctionReturn(0);

    PetscErrorCode ierr;
    PetscScalar dt, dt_scal, dt_fs, time_fs, Exx, Eyy, Rxx, Ryy;
    PetscInt ind, i, j, tnodes, tnodes_ori, L, sx, sy, sz, nx, ny, nz, step_fs, nx_solve, ny_solve;
    PetscScalar *topo_fs = PETSC_NULL, *topo_solve = PETSC_NULL, *topo_solve_refine = PETSC_NULL, *topo_alloc = PETSC_NULL; 
    PetscScalar *vx_alloc = PETSC_NULL, *vy_alloc = PETSC_NULL, *vz_alloc = PETSC_NULL;
    PetscScalar *vx_solve = PETSC_NULL, *vy_solve = PETSC_NULL, *vz_solve = PETSC_NULL;
    PetscScalar *vx_solve_refine = PETSC_NULL, *vy_solve_refine = PETSC_NULL, *vz_solve_refine = PETSC_NULL; 
    PetscScalar ***topo;

    // load global nx, ny, dt, time, rangeX, rangeY
    // nx, ny, dt, time
    FDSTAG   *fs;
    TSSol    *ts;
    Scaling  *scal;
    JacRes   *jr;   

    jr      = surf->jr;
    fs      = jr->fs;
    ts      = jr->ts;
    scal    = ts->scal;

    // load time
    dt      = ts->dt;
    dt_scal = dt * scal->time;
    dt_fs   = dt_scal * scal->time_fs; // (Myr) in LaMEM to (yr) in FastScape (GEO)
    time_fs = ts->time * scal->time + dt_scal; // time after finishing surface processes
    step_fs = ts->istep;

    // Gather topography and velocity
    tnodes_ori    = FSLib->nx_fs * FSLib->ny_fs;  

    topo_fs       = (PetscScalar *)malloc(tnodes_ori * sizeof(PetscScalar)); 
    topo_alloc    = (PetscScalar *)malloc(tnodes_ori * sizeof(PetscScalar));
    vx_alloc      = (PetscScalar *)malloc(tnodes_ori * sizeof(PetscScalar));  
    vy_alloc      = (PetscScalar *)malloc(tnodes_ori * sizeof(PetscScalar));  
    vz_alloc      = (PetscScalar *)malloc(tnodes_ori * sizeof(PetscScalar));    

    ierr = GatherVariableFromLaMEM(FSLib, topo_alloc, vx_alloc, vy_alloc, vz_alloc, step_fs);  CHKERRQ(ierr);        

    if(ISRankZero(PETSC_COMM_WORLD))
    {
        PetscScalar dt_max = FSLib->Max_dt; // Maximum step length, if dt_LaMEM is larger than this, use this
        PetscScalar dt_n = 0; //dt_residual
        PetscScalar quotient = dt_fs/dt_max;
        PetscInt nsteps = floor( dt_fs/dt_max );
        PetscScalar *topo_pass = PETSC_NULL, *vx_pass = PETSC_NULL, *vy_pass = PETSC_NULL, *vz_pass = PETSC_NULL; 

        // store the phase that is being sedimented
        surf->phase = FSLib->sedPhases;

        tnodes      = FSLib->nodes_solve;  
        nx_solve    = FSLib->nx_solve;
        ny_solve    = FSLib->ny_solve;

        topo_pass          = (PetscScalar *)malloc(tnodes * sizeof(PetscScalar));
        vx_pass            = (PetscScalar *)malloc(tnodes * sizeof(PetscScalar)); 
        vy_pass            = (PetscScalar *)malloc(tnodes * sizeof(PetscScalar));         
        vz_pass            = (PetscScalar *)malloc(tnodes * sizeof(PetscScalar));             

        // timestep
        if(1 > nsteps) 
        {
            nsteps = 1;
            dt_max = dt_fs;
        }
        else if( (1 <= nsteps) && (quotient != nsteps) )
        {
            nsteps = nsteps + 1;
            dt_n   = dt_fs - (nsteps - 1) * dt_max;
        }

        // get background strain rates
        ierr = BCGetBGStrainRates(jr->bc, &Exx, &Eyy, NULL, NULL, NULL, NULL, &Rxx, &Ryy, NULL);  CHKERRQ(ierr);

        if( Exx || Eyy )
        {
            ierr = FastScapeStretchGrid(FSLib);       CHKERRQ(ierr);  
        }

        // run fastscape when using 3D geodynamic model
        if(0 == FSLib->fs2D)
        {

            // for non uniform grid
            if( 1 == FSLib->non_uniform_grid )
            {
                // topography
                if( 0 == step_fs)
                {
                    ierr = InterpolationFor3DNonUniformGrid(FSLib, topo_alloc, 1);              CHKERRQ(ierr);  
                }

                // vx
                ierr = InterpolationFor3DNonUniformGrid(FSLib, vx_alloc, 2);                    CHKERRQ(ierr);             

                // vy
                ierr = InterpolationFor3DNonUniformGrid(FSLib, vy_alloc, 3);                    CHKERRQ(ierr);  
            
                // vz
                ierr = InterpolationFor3DNonUniformGrid(FSLib, vz_alloc, 4);                    CHKERRQ(ierr);  
            }               

            // don't apply refinement
            if( 1 == FSLib->refine)
            {
                if(0 == FSLib->non_uniform_grid)
                {
                    ierr = VecGetArray(FSLib->gtopo_fs, &topo_fs);                        CHKERRQ(ierr);

                    for(j = 0; j < ny_solve; j++) 
                    {   
                        for(i = 0; i < nx_solve; i++) 
                        {
                            ind            = j * nx_solve + i;


                            if (0 == step_fs)
                            {
                                topo_pass[ind] = topo_alloc[ind] * scal->length * scal->length_fs; // (km) in LaMEM to (m) in FastScape (GEO)
                            }
                            else
                            {
                                topo_pass[ind] = topo_fs[ind] * scal->length_fs;
                            }
                            vx_pass[ind]   = vx_alloc[ind] * scal->velocity / scal->velocity_fs; // (cm/yr) in LaMEM to (m/yr) in FastScape (GEO)
                            vy_pass[ind]   = vy_alloc[ind] * scal->velocity / scal->velocity_fs; // (cm/yr) in LaMEM to (m/yr) in FastScape (GEO)                            
                            vz_pass[ind]   = vz_alloc[ind] * scal->velocity / scal->velocity_fs; // (cm/yr) in LaMEM to (m/yr) in FastScape (GEO)
                        } 
                    }

                    ierr = VecRestoreArray(FSLib->gtopo_fs, &topo_fs);                        CHKERRQ(ierr);                     
                    
                    // run FastScape 
                    ierr = FastScapeFortranCppAdvc(FSLib, dt_max, dt_n, nsteps, step_fs, vx_pass, vy_pass, vz_pass, topo_pass); CHKERRQ(ierr);      
                }
                else
                {
                    ierr = VecGetArray(FSLib->gtopo_nug, &topo_solve);  CHKERRQ(ierr);
                    ierr = VecGetArray(FSLib->vx_nug, &vx_solve);       CHKERRQ(ierr);
                    ierr = VecGetArray(FSLib->vy_nug, &vy_solve);       CHKERRQ(ierr);                    
                    ierr = VecGetArray(FSLib->vz_nug, &vz_solve);       CHKERRQ(ierr);
                    
                    for(j = 0; j < ny_solve; j++) 
                    {   
                        for(i = 0; i < nx_solve; i++) 
                        {
                            ind            = j * nx_solve + i;

                            topo_pass[ind] = topo_solve[ind] * scal->length_fs; // (km) in LaMEM to (m) in FastScape (GEO)
                            vx_pass[ind]   = vx_solve[ind] * scal->velocity / scal->velocity_fs; // (cm/yr) in LaMEM to (m/yr) in FastScape (GEO)
                            vy_pass[ind]   = vy_solve[ind] * scal->velocity / scal->velocity_fs; // (cm/yr) in LaMEM to (m/yr) in FastScape (GEO)  
                            vz_pass[ind]   = vz_solve[ind] * scal->velocity / scal->velocity_fs; // (cm/yr) in LaMEM to (m/yr) in FastScape (GEO)
                        } 
                    }

                    ierr = VecRestoreArray(FSLib->gtopo_nug, &topo_solve);  CHKERRQ(ierr);
                    ierr = VecRestoreArray(FSLib->vx_nug, &vx_solve);       CHKERRQ(ierr);
                    ierr = VecRestoreArray(FSLib->vy_nug, &vy_solve);       CHKERRQ(ierr);                    
                    ierr = VecRestoreArray(FSLib->vz_nug, &vz_solve);       CHKERRQ(ierr);
    
                    // run FastScape 
                    ierr = FastScapeFortranCppAdvc(FSLib, dt_max, dt_n, nsteps, step_fs, vx_pass, vy_pass, vz_pass, topo_pass); CHKERRQ(ierr);        
                }                      
            }                           
            // apply refinement
            else 
            {
                printf("Refined times                    : %d\n", FSLib->refine);
                printf("Refined grid cells [nx, ny]      : [%d, %d] \n", FSLib->nx_refine, FSLib->ny_refine);

                ierr = VecGetArray(FSLib->gtopo_refine, &topo_solve_refine);  CHKERRQ(ierr);
                ierr = VecGetArray(FSLib->vx_refine, &vx_solve_refine);       CHKERRQ(ierr);
                ierr = VecGetArray(FSLib->vy_refine, &vy_solve_refine);       CHKERRQ(ierr);
                ierr = VecGetArray(FSLib->vz_refine, &vz_solve_refine);       CHKERRQ(ierr);
    
                if(0 == FSLib->non_uniform_grid)
                {
                    // load velocity and topography field
                    if(0 == step_fs)
                    {
                        ierr = BilinearInterpolate(FSLib, topo_alloc, topo_solve_refine, topo_pass, scal, 1, FSLib->nx_refine, FSLib->ny_refine); CHKERRQ(ierr);
                    }
                    else
                    {
                        for(j = 0; j < ny_solve; j++)
                        {
                            for(i = 0; i < nx_solve; i++)
                            {
                                ind = j * nx_solve + i;

                                topo_pass[ind] = topo_solve_refine[ind] * scal->length_fs;
                            }
                        }
                    }
                    ierr = BilinearInterpolate(FSLib, vx_alloc, vx_solve_refine, vx_pass, scal, 2, FSLib->nx_refine, FSLib->ny_refine); CHKERRQ(ierr);
                    ierr = BilinearInterpolate(FSLib, vy_alloc, vy_solve_refine, vy_pass, scal, 2, FSLib->nx_refine, FSLib->ny_refine); CHKERRQ(ierr);                    
                    ierr = BilinearInterpolate(FSLib, vz_alloc, vz_solve_refine, vz_pass, scal, 2, FSLib->nx_refine, FSLib->ny_refine); CHKERRQ(ierr);
    
                    // run FastScape
                    ierr = FastScapeFortranCppAdvc(FSLib, dt_max, dt_n, nsteps, step_fs, vx_pass, vy_pass, vz_pass, topo_pass); CHKERRQ(ierr);
                }
                else
                {
                    ierr = VecGetArray(FSLib->gtopo_nug, &topo_solve);  CHKERRQ(ierr);
                    ierr = VecGetArray(FSLib->vz_nug, &vz_solve);       CHKERRQ(ierr);
                    ierr = VecGetArray(FSLib->vx_nug, &vx_solve);       CHKERRQ(ierr);
                    ierr = VecGetArray(FSLib->vy_nug, &vy_solve);       CHKERRQ(ierr);
                    
                    // load velocity and topography field
                    if(0 == step_fs)
                    {
                        ierr = BilinearInterpolate(FSLib, topo_solve, topo_solve_refine, topo_pass, scal, 1, FSLib->nx_refine, FSLib->ny_refine); CHKERRQ(ierr);
                    }
                    else
                    {
                        for(j = 0; j < ny_solve; j++)
                        {
                            for(i = 0; i < nx_solve; i++)
                            {
                                ind = j * nx_solve + i;

                                topo_pass[ind] = topo_solve_refine[ind] * scal->length_fs;
                            }
                        }
                    }

                    ierr = BilinearInterpolate(FSLib, vx_solve, vx_solve_refine, vx_pass, scal, 2, FSLib->nx_refine, FSLib->ny_refine); CHKERRQ(ierr);
                    ierr = BilinearInterpolate(FSLib, vy_solve, vy_solve_refine, vy_pass, scal, 2, FSLib->nx_refine, FSLib->ny_refine); CHKERRQ(ierr);                    
                    ierr = BilinearInterpolate(FSLib, vz_solve, vz_solve_refine, vz_pass, scal, 2, FSLib->nx_refine, FSLib->ny_refine); CHKERRQ(ierr);                    
                }

                // run FastScape
                ierr = FastScapeFortranCppAdvc(FSLib, dt_max, dt_n, nsteps, step_fs, vx_pass, vy_pass, vz_pass, topo_pass); CHKERRQ(ierr);
                    
                ierr = VecRestoreArray(FSLib->gtopo_nug, &topo_solve);  CHKERRQ(ierr);
                ierr = VecRestoreArray(FSLib->vz_nug, &vz_solve);       CHKERRQ(ierr);
                ierr = VecRestoreArray(FSLib->vx_nug, &vx_solve);       CHKERRQ(ierr);
                ierr = VecRestoreArray(FSLib->vy_nug, &vy_solve);       CHKERRQ(ierr); 
                
                ierr = VecRestoreArray(FSLib->gtopo_refine, &topo_solve_refine);  CHKERRQ(ierr);
                ierr = VecRestoreArray(FSLib->vx_refine, &vx_solve_refine);       CHKERRQ(ierr);
                ierr = VecRestoreArray(FSLib->vy_refine, &vy_solve_refine);       CHKERRQ(ierr);
                ierr = VecRestoreArray(FSLib->vz_refine, &vz_solve_refine);       CHKERRQ(ierr);
            }
        }
        // run fastscape when using 2D geodynamic model
        else
        {
            ierr = VecGetArray(FSLib->gtopo_extend, &topo_solve);  CHKERRQ(ierr);
            ierr = VecGetArray(FSLib->vx_extend, &vx_solve);       CHKERRQ(ierr);
            ierr = VecGetArray(FSLib->vy_extend, &vy_solve);       CHKERRQ(ierr); 
            ierr = VecGetArray(FSLib->vz_extend, &vz_solve);       CHKERRQ(ierr); 

            if(0 == FSLib->non_uniform_grid)
            {
                // don't apply refinement
                if(1 == FSLib->refine)
                {
                    // load velocity and topography field
                    // step == 0, using advection in LaMEM, after that, using advection in FastScape (due to the initial value from FastScape but not from LaMEM)
                    // topography, advect in FastScape, Advect1D: advect in x or y direction, FastScape: advect in z direction
                    if( 0 == step_fs) 
                    {
                        ierr = Extended2D(FSLib, topo_alloc, topo_solve, topo_pass, scal, 1); CHKERRQ(ierr);
                    }
                    else
                    {
                        for(j = 0; j < ny_solve; j++)
                        {
                            for(i = 0; i < nx_solve; i++)
                            {
                                ind = j * nx_solve + i;

                                topo_pass[ind] = topo_solve[ind] * scal->length_fs;
                            }
                        }
                    }

                    // velocity
                    // vx, vy
                    if(1 == FSLib->extendedX)
                    {
                        ierr = Extended2D(FSLib, vy_alloc, vy_solve, vy_pass, scal, 2); CHKERRQ(ierr);

                        for(j = 0; j < ny_solve; j++)
                        {
                            for(i = 0; i < nx_solve; i++)
                            {
                                ind = j * nx_solve + i;

                                vx_pass[ind]  = 0;
                            }
                        }
                    }
                    else
                    {
                        ierr = Extended2D(FSLib, vx_alloc, vx_solve, vx_pass, scal, 2); CHKERRQ(ierr);    

                        for(i = 0; i < nx_solve; i++)
                        {
                            for(j = 0; j < ny_solve; j++)
                            {
                                ind = j * nx_solve + i;

                                vy_pass[ind]  = 0;
                            }
                        }
                    }
                    // vz
                    ierr = Extended2D(FSLib, vz_alloc, vz_solve, vz_pass, scal, 2); CHKERRQ(ierr);

                    // run FastScape
                    ierr = FastScapeFortranCppAdvc(FSLib, dt_max, dt_n, nsteps, step_fs, vx_pass, vy_pass, vz_pass, topo_pass); CHKERRQ(ierr);
                }
                // apply refinement
                else
                {
                    // load velocity and topography field
                    ierr = VecGetArray(FSLib->gtopo_et_refine, &topo_solve_refine);   CHKERRQ(ierr);
                    ierr = VecGetArray(FSLib->vx_et_refine, &vx_solve_refine);        CHKERRQ(ierr);
                    ierr = VecGetArray(FSLib->vy_et_refine, &vy_solve_refine);        CHKERRQ(ierr);
                    ierr = VecGetArray(FSLib->vz_et_refine, &vz_solve_refine);        CHKERRQ(ierr);

                    if( 0 == step_fs) 
                    {
                        ierr = Extended2DRefine(FSLib, topo_alloc, topo_solve);   CHKERRQ(ierr);
                        ierr = BilinearInterpolate(FSLib, topo_solve, topo_solve_refine, topo_pass, scal, 1, FSLib->etRefineXNodes, FSLib->etRefineYNodes); CHKERRQ(ierr);
                    }
                    else
                    {
                        for(j = 0; j < ny_solve; j++)
                        {
                            for(i = 0; i < nx_solve; i++)
                            {
                                ind = j * nx_solve + i;

                                topo_pass[ind] = topo_solve_refine[ind] * scal->length_fs;
                            }
                        }
                    }

                    if(1 == FSLib->extendedX)
                    {
                        ierr = Extended2DRefine(FSLib, vy_alloc, vy_solve);           CHKERRQ(ierr);

                        for(j = 0; j < ny_solve; j++)
                        {
                            for(i = 0; i < nx_solve; i++)
                            {
                                ind = j * nx_solve + i;

                                vx_solve[ind]  = 0;
                            }
                        }
                    }
                    else
                    {
                        ierr = Extended2DRefine(FSLib, vx_alloc, vx_solve);           CHKERRQ(ierr);   

                        for(i = 0; i < nx_solve; i++)
                        {
                            for(j = 0; j < ny_solve; j++)
                            {
                                ind = j * nx_solve + i;

                                vy_solve[ind]  = 0;
                            }
                        }
                    }
                        
                    ierr = Extended2DRefine(FSLib, vz_alloc, vz_solve);              CHKERRQ(ierr);

                    ierr = BilinearInterpolate(FSLib, vx_solve, vx_solve_refine, vx_pass, scal, 2, FSLib->etRefineXNodes, FSLib->etRefineYNodes); CHKERRQ(ierr);
                    ierr = BilinearInterpolate(FSLib, vy_solve, vy_solve_refine, vy_pass, scal, 2, FSLib->etRefineXNodes, FSLib->etRefineYNodes); CHKERRQ(ierr);                    
                    ierr = BilinearInterpolate(FSLib, vz_solve, vz_solve_refine, vz_pass, scal, 2, FSLib->etRefineXNodes, FSLib->etRefineYNodes); CHKERRQ(ierr);

                    ierr = VecRestoreArray(FSLib->gtopo_et_refine, &topo_solve_refine);  CHKERRQ(ierr);
                    ierr = VecRestoreArray(FSLib->vx_et_refine, &vx_solve_refine);       CHKERRQ(ierr);
                    ierr = VecRestoreArray(FSLib->vy_et_refine, &vy_solve_refine);       CHKERRQ(ierr);
                    ierr = VecRestoreArray(FSLib->vz_et_refine, &vz_solve_refine);       CHKERRQ(ierr);            

                    // run FastScape
                    ierr = FastScapeFortranCppAdvc(FSLib, dt_max, dt_n, nsteps, step_fs, vx_pass, vy_pass, vz_pass, topo_pass); CHKERRQ(ierr);
                }
            }
            else
            {
                // don't apply refinement
                if(1 == FSLib->refine)
                {
                    // load velocity and topography field
                    // step == 0, using advection in LaMEM, after that, using advection in FastScape (due to the initial value from FastScape but not from LaMEM)
                    // topography, advect in FastScape, Advect1D: advect in x or y direction, FastScape: advect in z direction
                    if( 0 == step_fs) 
                    {
                        ierr = Extended2D(FSLib, topo_alloc, topo_solve, topo_pass, scal, 1); CHKERRQ(ierr);
                    }
                    else
                    {
                        for(j = 0; j < ny_solve; j++)
                        {
                            for(i = 0; i < nx_solve; i++)
                            {
                                ind = j * nx_solve + i;

                                topo_pass[ind] = topo_solve[ind] * scal->length_fs;
                            }
                        }
                    }

                    // velocity
                    // vx, vy
                    if(1 == FSLib->extendedX)
                    {
                        ierr = Extended2D(FSLib, vy_alloc, vy_solve, vy_pass, scal, 2); CHKERRQ(ierr);

                        for(j = 0; j < ny_solve; j++)
                        {
                            for(i = 0; i < nx_solve; i++)
                            {
                                ind = j * nx_solve + i;

                                vx_pass[ind]  = 0;
                            }
                        }
                    }
                    else
                    {
                        ierr = Extended2D(FSLib, vx_alloc, vx_solve, vx_pass, scal, 2); CHKERRQ(ierr);    

                        for(i = 0; i < nx_solve; i++)
                        {
                            for(j = 0; j < ny_solve; j++)
                            {
                                ind = j * nx_solve + i;

                                vy_pass[ind]  = 0;
                            }
                        }
                    }
                    // vz
                    ierr = Extended2D(FSLib, vz_alloc, vz_solve, vz_pass, scal, 2); CHKERRQ(ierr);

                    // run FastScape
                    ierr = FastScapeFortranCppAdvc(FSLib, dt_max, dt_n, nsteps, step_fs, vx_pass, vy_pass, vz_pass, topo_pass); CHKERRQ(ierr);
                }
                // apply refinement
                else
                {
                    // load velocity and topography field
                    ierr = VecGetArray(FSLib->gtopo_et_refine, &topo_solve_refine);   CHKERRQ(ierr);
                    ierr = VecGetArray(FSLib->vx_et_refine, &vx_solve_refine);        CHKERRQ(ierr);
                    ierr = VecGetArray(FSLib->vy_et_refine, &vy_solve_refine);        CHKERRQ(ierr);
                    ierr = VecGetArray(FSLib->vz_et_refine, &vz_solve_refine);        CHKERRQ(ierr);

                    if( 0 == step_fs) 
                    {
                        ierr = Extended2DRefine(FSLib, topo_alloc, topo_solve);   CHKERRQ(ierr);
                        ierr = BilinearInterpolate(FSLib, topo_solve, topo_solve_refine, topo_pass, scal, 1, FSLib->etRefineXNodes, FSLib->etRefineYNodes); CHKERRQ(ierr);
                    }
                    else
                    {
                        for(j = 0; j < ny_solve; j++) 
                        {
                            for(i = 0; i < nx_solve; i++)
                            {
                                ind = j * nx_solve + i;

                                topo_pass[ind] = topo_solve_refine[ind] * scal->length_fs;
                            }
                        }
                    }

                    if(1 == FSLib->extendedX)
                    {
                        ierr = Extended2DRefine(FSLib, vy_alloc, vy_solve);           CHKERRQ(ierr);

                        for(j = 0; j < ny_solve; j++)
                        {
                            for(i = 0; i < nx_solve; i++)
                            {
                                ind = j * nx_solve + i;

                                vx_solve[ind]  = 0;
                            }
                        }
                    }
                    else
                    {
                        ierr = Extended2DRefine(FSLib, vx_alloc, vx_solve);           CHKERRQ(ierr);   

                        for(i = 0; i < nx_solve; i++)
                        {
                            for(j = 0; j < ny_solve; j++)
                            {
                                ind = j * nx_solve + i;

                                vy_solve[ind]  = 0;
                            }
                        }
                    }
                        
                    ierr = Extended2DRefine(FSLib, vz_alloc, vz_solve);              CHKERRQ(ierr);

                    ierr = BilinearInterpolate(FSLib, vx_solve, vx_solve_refine, vx_pass, scal, 2, FSLib->etRefineXNodes, FSLib->etRefineYNodes); CHKERRQ(ierr);
                    ierr = BilinearInterpolate(FSLib, vy_solve, vy_solve_refine, vy_pass, scal, 2, FSLib->etRefineXNodes, FSLib->etRefineYNodes); CHKERRQ(ierr);                    
                    ierr = BilinearInterpolate(FSLib, vz_solve, vz_solve_refine, vz_pass, scal, 2, FSLib->etRefineXNodes, FSLib->etRefineYNodes); CHKERRQ(ierr);

                    ierr = VecRestoreArray(FSLib->gtopo_et_refine, &topo_solve_refine);  CHKERRQ(ierr);
                    ierr = VecRestoreArray(FSLib->vx_et_refine, &vx_solve_refine);       CHKERRQ(ierr);
                    ierr = VecRestoreArray(FSLib->vy_et_refine, &vy_solve_refine);       CHKERRQ(ierr);
                    ierr = VecRestoreArray(FSLib->vz_et_refine, &vz_solve_refine);       CHKERRQ(ierr);
                    
                    // run FastScape
                    ierr = FastScapeFortranCppAdvc(FSLib, dt_max, dt_n, nsteps, step_fs, vx_pass, vy_pass, vz_pass, topo_pass); CHKERRQ(ierr);
                }
            }

            ierr = VecRestoreArray(FSLib->gtopo_extend, &topo_solve);  CHKERRQ(ierr);
            ierr = VecRestoreArray(FSLib->vx_extend, &vx_solve);       CHKERRQ(ierr);
            ierr = VecRestoreArray(FSLib->vy_extend, &vy_solve);       CHKERRQ(ierr);
            ierr = VecRestoreArray(FSLib->vz_extend, &vz_solve);       CHKERRQ(ierr);       
        }
        
        ierr = VecGetArray(FSLib->gtopo_fs, &topo_fs);                        CHKERRQ(ierr);

        // pass data to original LaMEM grid
        if(0 == FSLib->fs2D)
        {
            ierr = PassValue3D(FSLib, topo_pass, topo_fs);                CHKERRQ(ierr); 

        }
        else
        {
            if(1 == FSLib->refine)
            {
                ierr = PassValue2D(FSLib, topo_pass, topo_fs);         CHKERRQ(ierr); 
            }
            else
            {
                ierr = Pass2DValueRefine(FSLib, topo_pass, topo_fs); CHKERRQ(ierr); 
            }
        }

        // save data in a new grid used by fastscape
        if(0 == step_fs || 0 == (step_fs + 1) % FSLib->surf_out_nstep)
        {
            ierr = FastScapeSave(FSLib, step_fs, time_fs);              CHKERRQ(ierr);
        }

        free(topo_pass);
        free(vx_pass); 
        free(vy_pass); 
        free(vz_pass); 

        topo_pass          = PETSC_NULL;
        vx_pass            = PETSC_NULL;
        vy_pass            = PETSC_NULL;        
        vz_pass            = PETSC_NULL;
    }
       
    // Broadcast      
    if(ISParallel(PETSC_COMM_WORLD))
    {
        ierr = MPI_Bcast(topo_fs, (PetscMPIInt)tnodes_ori, MPIU_SCALAR, (PetscMPIInt)0, PETSC_COMM_WORLD);   CHKERRQ(ierr);    
    }

    L    = (PetscInt)fs->dsz.rank;
    ierr = DMDAGetCorners(fs->DA_COR, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(surf->DA_SURF, surf->gtopo, &topo);       CHKERRQ(ierr);

    // Save topography in different rank
    START_PLANE_LOOP
    {
        topo[L][j][i] = topo_fs[j * FSLib->nx_fs + i] / scal->length; // (km)
    }
    END_PLANE_LOOP  

    ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->gtopo, &topo);  CHKERRQ(ierr);
    ierr = VecRestoreArray(FSLib->gtopo_fs, &topo_fs);              CHKERRQ(ierr);

    // compute ghosted version of the topography
    GLOBAL_TO_LOCAL(surf->DA_SURF, surf->gtopo, surf->ltopo);

    // compute & store average topography
    ierr = FreeSurfGetAvgTopo(surf);                                CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");
    
    free(topo_alloc);
    free(vz_alloc);
    free(vx_alloc);
    free(vy_alloc);
    free(topo_fs);

    topo_alloc = PETSC_NULL;
    vz_alloc = PETSC_NULL;
    vx_alloc = PETSC_NULL;
    vy_alloc = PETSC_NULL;
    topo_fs = PETSC_NULL;

    PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PassValue2D(FastScapeLib *FSLib, PetscScalar *topo_pass_f, PetscScalar *topo_fs)
{
    PetscErrorCode ierr;
    PetscScalar sum_topo;
    PetscInt i, j, ind; 
    PetscScalar *gtopo_extend = PETSC_NULL;
    PetscScalar *topo_aver = PETSC_NULL;
    PetscScalar *topo_aver_ori = PETSC_NULL;
    FSGrid  *fsX;
    FSGrid  *fsY;
    Scaling *scal;
    PetscInt m, n, mm, nn, ind_a, ind_b;
    PetscScalar dx, dy, wtx, wty, x_coor, y_coor, x_begin = 0.0, y_begin = 0.0;

    fsX     = &FSLib->fsX;
    fsY     = &FSLib->fsY;
    scal    = FSLib->scal;

    if(1 == FSLib->non_uniform_grid)
    {    
        if(1 == FSLib->extendedX)
        {
            y_begin = FSLib->ncoor_ori_y[0];
        }
        else
        {
            x_begin = FSLib->ncoor_ori_x[0];
        }
    }

    ierr = VecGetArray(FSLib->gtopo_extend,  &gtopo_extend);       CHKERRQ(ierr);
        
    // extend in x-direction
    if(1 == FSLib->extendedX)
    {
        topo_aver     = (PetscScalar *)malloc(FSLib->extendedYNodes * sizeof(PetscScalar));
        topo_aver_ori = (PetscScalar *)malloc(FSLib->ny_fs          * sizeof(PetscScalar));

        for(j = 0; j < FSLib->extendedYNodes; j++)
        {
            sum_topo = 0;

            for(i = 0; i < FSLib->extendedXNodes; i++)
            {
                ind = j * FSLib->extendedXNodes + i;

                // correction
                gtopo_extend[ind] = topo_pass_f[ind] / scal->length_fs;        // m -> km        

                if(gtopo_extend[ind] > FSLib->rangeZ_end)   gtopo_extend[ind] = FSLib->rangeZ_end;
                if(gtopo_extend[ind] < FSLib->rangeZ_begin) gtopo_extend[ind] = FSLib->rangeZ_begin;

                sum_topo          = sum_topo + gtopo_extend[ind]; 
            }

            topo_aver[j] = sum_topo / FSLib->extendedXNodes;
        }

        if(0 == FSLib->non_uniform_grid)
        {   
            for(j = 0; j < FSLib->ny_fs; j++)
            {
                for(i = 0; i < FSLib->nx_fs; i++)
                {
                    ind          = j * FSLib->nx_fs + i;

                    topo_fs[ind] = topo_aver[j];         
                }         
            }
        }
        else
        {
            for(j = 0; j < FSLib->ny_fs; j++)
            {  
                y_coor = FSLib->ncoor_ori_y[j];
                dy     = fsY->dx;

                // get nearest 2 index
                // y-direction
                n  = floor( (y_coor - y_begin) / fsY->dx );
                nn = n + 1;
                if( FSLib->extendedYNodes == nn) nn -= 1;

                // interpolate
                ind = j;    
                ind_a = n;
                ind_b = nn; 

                // linear interpolation
                wty   = y_coor - fsY->ncoor_extend[n];
                 
                if(n == nn)
                {
                    topo_aver_ori[ind] = topo_aver[ind_a];
                }
                else
                {
                    topo_aver_ori[ind] =   ( (dy - wty) / dy ) * topo_aver[ind_a]+ ( wty / dy ) * topo_aver[ind_b];        
                }
            }  

            for(j = 0; j < FSLib->ny_fs; j++)
            {
                for(i = 0; i < FSLib->nx_fs; i++)
                {
                    ind          = j * FSLib->nx_fs + i;

                    topo_fs[ind] = topo_aver_ori[j];         
                }         
            }            
        }
    }
    // extend in y-direction
    else
    {
        topo_aver     = (PetscScalar *)malloc(FSLib->extendedXNodes  * sizeof(PetscScalar));
        topo_aver_ori = (PetscScalar *)malloc(FSLib->nx_fs           * sizeof(PetscScalar));

        for(i = 0; i < FSLib->extendedXNodes; i++)
        {
            sum_topo = 0;

            for(j = 0; j < FSLib->extendedYNodes; j++)
            {
                ind = j * FSLib->extendedXNodes + i;

                gtopo_extend[ind] = topo_pass_f[ind] / scal->length_fs;

                // correction
                if(gtopo_extend[ind] > FSLib->rangeZ_end )   gtopo_extend[ind] = FSLib->rangeZ_end;
                if(gtopo_extend[ind] < FSLib->rangeZ_begin ) gtopo_extend[ind] = FSLib->rangeZ_begin;

                sum_topo          = sum_topo + gtopo_extend[ind];
            }

            topo_aver[i] = sum_topo/FSLib->extendedYNodes;         
        }

        if(0 == FSLib->non_uniform_grid)
        {   
            for(i = 0; i < FSLib->nx_fs; i++)
            {
                for(j = 0; j < FSLib->ny_fs; j++)
                {
                    ind          = j * FSLib->nx_fs + i;
    
                    topo_fs[ind] = topo_aver[i];  
                }         
            }  
        }
        else
        {
            for(i = 0; i < FSLib->nx_fs; i++)
            {  
                x_coor = FSLib->ncoor_ori_x[i];
                dx     = fsX->dx;

                // get nearest 2 index
                // x-direction
                m  = floor( (x_coor - x_begin) / dx );
                mm = m + 1;
                if( FSLib->extendedXNodes == mm) mm -= 1;

                // interpolate
                ind = i;    
                ind_a = m;
                ind_b = mm; 

                // linear interpolation
                wtx   = x_coor - fsX->ncoor_extend[m];

                if(m == mm)
                {
                    topo_aver_ori[ind] = topo_aver[ind_a];
                }
                else           
                {                                                     
                    topo_aver_ori[ind] =   ( (dx - wtx) / dx ) * topo_aver[ind_a]+ ( wtx / dx ) * topo_aver[ind_b];      
                }        
            }  

            for(i = 0; i < FSLib->nx_fs; i++)
            {
                for(j = 0; j < FSLib->ny_fs; j++)
                {
                    ind          = j * FSLib->nx_fs + i;
    
                    topo_fs[ind] = topo_aver_ori[i];  
                }         
            }             
        }
    }

    ierr = VecRestoreArray(FSLib->gtopo_extend, &gtopo_extend);  CHKERRQ(ierr);   

    free(topo_aver);
    free(topo_aver_ori);
    topo_aver = PETSC_NULL;
    topo_aver_ori = PETSC_NULL;

    PetscFunctionReturn(0);    
}
//---------------------------------------------------------------------------
PetscErrorCode Pass2DValueRefine(FastScapeLib *FSLib, PetscScalar *topo_pass_f, PetscScalar *topo_fs)
{
    PetscErrorCode ierr;
    PetscScalar sum_topo;
    PetscInt i, j, ind, ind2, count, countX = 0, countY = 0;

    PetscScalar *topo_et_refine = PETSC_NULL, *topo_extend = PETSC_NULL;
    PetscScalar *topo_aver = PETSC_NULL, *topo_aver_ori = PETSC_NULL;
    FSGrid  *fsX;
    FSGrid  *fsY;
    Scaling *scal;
    PetscInt m, n, mm, nn, ind_a, ind_b;
    PetscScalar dx, dy, wtx, wty, x_coor, y_coor, x_begin = 0.0, y_begin = 0.0;

    fsX     = &FSLib->fsX;
    fsY     = &FSLib->fsY;
    scal    = FSLib->scal;

    if(1 == FSLib->non_uniform_grid)
    {
        if(1 == FSLib->extendedX)
        {
            y_begin = FSLib->ncoor_ori_y[0];
        }
        else
        {
            x_begin = FSLib->ncoor_ori_x[0];
        }
    }
    
    ierr = VecGetArray(FSLib->gtopo_et_refine, &topo_et_refine);  CHKERRQ(ierr);
    ierr = VecGetArray(FSLib->gtopo_extend, &topo_extend);        CHKERRQ(ierr); 

    for(j = 0; j < FSLib->etRefineYNodes; j++)
    {
        for(i = 0; i < FSLib->etRefineXNodes; i++)
        {
            ind = j * FSLib->etRefineXNodes + i;

            topo_et_refine[ind] = topo_pass_f[ind] / scal->length_fs;

            // correction
            if(topo_et_refine[ind] > FSLib->rangeZ_end)   topo_et_refine[ind] = FSLib->rangeZ_end;
            if(topo_et_refine[ind] < FSLib->rangeZ_begin) topo_et_refine[ind] = FSLib->rangeZ_begin;
        }
    }

    // extend in x-direction
    if(1 == FSLib->extendedX)
    {
        topo_aver_ori = (PetscScalar *)malloc(FSLib->ny_fs * sizeof(PetscScalar));
        topo_aver = (PetscScalar *)malloc(FSLib->extendedYNodes * sizeof(PetscScalar));

        count = 0; 

        for(j = 0; j < FSLib->etRefineYNodes; j += FSLib->refine) 
        {  
            sum_topo = 0;

            for(i = 0; i < FSLib->etRefineXNodes; i += FSLib->refine) 
            {
                ind = j * FSLib->etRefineXNodes+i; 

                if(0 == ind % FSLib->etRefineXNodes) 
                { 
                    countX = 0;
                    if(0 != ind) 
                    {
                        countY++; 
                    }
                }

                ind2              = countY * FSLib->extendedXNodes + countX;                

                sum_topo          = sum_topo + topo_et_refine[ind]; // m -> km
                topo_extend[ind2] = topo_et_refine[ind];

                countX++;

            }

            topo_aver[count] = sum_topo / FSLib->extendedXNodes; 

            count++;           
        }

        if(0 == FSLib->non_uniform_grid)
        {
            for(j = 0; j < FSLib->ny_fs; j++)
            {
                for(i = 0; i < FSLib->nx_fs; i++)
                {
                    ind          = j * FSLib->nx_fs + i;
                    topo_fs[ind] = topo_aver[j];             
                }         
            }        
        }
        else
        {
            for(j = 0; j < FSLib->ny_fs; j++)
            {  
                y_coor = FSLib->ncoor_ori_y[j];
                dy     = fsY->dx;

                // get nearest 2 index
                // y-direction
                n  = floor( (y_coor - y_begin) / dy );
                nn = n + 1;
                if( FSLib->extendedYNodes == nn) nn -= 1;

                // interpolate
                ind = j;    
                ind_a = n;
                ind_b = nn; 

                // linear interpolation
                wty   = y_coor - fsY->ncoor_extend[n];
                           
                if(n == nn)
                {
                    topo_aver_ori[ind] = topo_aver[ind_a];
                }
                else
                {
                    topo_aver_ori[ind] =   ( (dy - wty) / dy ) * topo_aver[ind_a]+ ( wty / dy ) * topo_aver[ind_b];   
                }     
            }  

            for(j = 0; j < FSLib->ny_fs; j++)
            {
                for(i = 0; i < FSLib->nx_fs; i++)
                {
                    ind          = j * FSLib->nx_fs + i;

                    topo_fs[ind] = topo_aver_ori[j];         
                }         
            }            
        }
    }
    // extend in y-direction
    else
    {
        topo_aver_ori = (PetscScalar *)malloc(FSLib->nx_fs          * sizeof(PetscScalar));
        topo_aver     = (PetscScalar *)malloc(FSLib->extendedXNodes * sizeof(PetscScalar));

        count = 0; 

        for(i = 0; i < FSLib->etRefineXNodes; i += FSLib->refine)
        {
            sum_topo = 0;

            for(j = 0; j < FSLib->etRefineYNodes; j += FSLib->refine)
            {
                ind = j * FSLib->etRefineXNodes+i; 

                if(0 == j % FSLib->etRefineYNodes) 
                {
                    countY = 0;

                    if(0 != ind) 
                    {
                        countX++; 
                    }
                }

                ind2              = countY * FSLib->extendedXNodes + countX;                

                sum_topo          = sum_topo + topo_et_refine[ind]; // m -> km
                topo_extend[ind2] = topo_et_refine[ind];

                countY++;
            }
            topo_aver[count] = sum_topo/FSLib->extendedYNodes;         

            count++;     
        }

        if(0 == FSLib->non_uniform_grid)
        {   
            for(i = 0; i < FSLib->nx_fs; i++)
            {
                for(j = 0; j < FSLib->ny_fs; j++)
                {
                    ind          = j * FSLib->nx_fs + i;

                    topo_fs[ind] = topo_aver[i];                    
                }         
            }  
        }
        else
        {
            for(i = 0; i < FSLib->nx_fs; i++)
            {  
                x_coor = FSLib->ncoor_ori_x[i];
                dx     = fsX->dx;

                // get nearest 2 index
                // x-direction
                m  = floor( (x_coor - x_begin) / dx );
                mm = m + 1;
                if( FSLib->extendedXNodes == mm) mm -= 1;

                // interpolate
                ind = i;    
                ind_a = m;
                ind_b = mm; 

                // linear interpolation
                wtx   = x_coor - fsX->ncoor_extend[m];
                                               
                if(m == mm)
                {
                    topo_aver_ori[ind] = topo_aver[ind_a];
                }
                else
                {
                    topo_aver_ori[ind] =   ( (dx - wtx) / dx ) * topo_aver[ind_a] + ( wtx / dx ) * topo_aver[ind_b];
                }              
            }  

            for(i = 0; i < FSLib->nx_fs; i++)
            {
                for(j = 0; j < FSLib->ny_fs; j++)
                {
                    ind          = j * FSLib->nx_fs + i;
    
                    topo_fs[ind] = topo_aver_ori[i];  
                }         
            }             
        }
    }

    ierr = VecRestoreArray(FSLib->gtopo_et_refine,  &topo_et_refine);  CHKERRQ(ierr);
    ierr = VecRestoreArray(FSLib->gtopo_extend,  &topo_extend);        CHKERRQ(ierr); 

    free(topo_aver);
    free(topo_aver_ori);
    topo_aver = PETSC_NULL;  
    topo_aver_ori = PETSC_NULL;


    PetscFunctionReturn(0);    
}
//---------------------------------------------------------------------------
PetscErrorCode PassValue3D(FastScapeLib *FSLib, PetscScalar *topo_pass_f, PetscScalar *topo_fs)
{
    PetscErrorCode ierr;
    PetscInt i, j, ind, ind2;
    PetscScalar *topo_refine = PETSC_NULL, *topo_nug = PETSC_NULL;
    
    FSGrid  *fsX;
    FSGrid  *fsY;
    Scaling *scal;

    PetscInt m, n, mm, nn, ind_a, ind_b, ind_c, ind_d;
    PetscScalar dx, dy, wtx, wty, x_coor, y_coor, x_begin = 0.0, y_begin = 0.0;


    fsX     = &FSLib->fsX;
    fsY     = &FSLib->fsY;
    scal    = FSLib->scal;
        
    if(1 == FSLib->non_uniform_grid)
    {
        x_begin = FSLib->ncoor_ori_x[0];
        y_begin = FSLib->ncoor_ori_y[0];
    }

    if(0 == FSLib->non_uniform_grid)
    {
        if(1 == FSLib->refine)
        {
            for(j = 0; j < FSLib->ny_fs; j++) 
            {   
                for(i = 0; i < FSLib->nx_fs; i++) 
                {
                    ind = j * FSLib->nx_fs + i;

                    topo_fs[ind] = topo_pass_f[ind] / scal->length_fs; // m(FastScape) to km(LaMEM) 
            
                    if(topo_fs[ind] > FSLib->rangeZ_end)   topo_fs[ind] = FSLib->rangeZ_end;
                    if(topo_fs[ind] < FSLib->rangeZ_begin) topo_fs[ind] = FSLib->rangeZ_begin;

                }
            }
        }            
        else
        {
            PetscInt countX = 0, countY = 0;

            ierr = VecGetArray(FSLib->gtopo_refine, &topo_refine);         CHKERRQ(ierr);

            for(j = 0; j < FSLib->ny_refine; j++)
            {
                for (i = 0; i < FSLib->nx_refine; i++)
                {
                    ind = j * FSLib->nx_refine + i;

                    topo_refine[ind] = topo_pass_f[ind] / scal->length_fs; // m(FastScape) to km(LaMEM)

                    if(topo_refine[ind] > FSLib->rangeZ_end)   topo_refine[ind] = FSLib->rangeZ_end;
                    if(topo_refine[ind] < FSLib->rangeZ_begin) topo_refine[ind] = FSLib->rangeZ_begin;
                }
            }

            for(j = 0; j < FSLib->ny_refine; j += FSLib->refine) 
            {   
                for(i = 0; i < FSLib->nx_refine; i += FSLib->refine) 
                {
                    ind = j * FSLib->nx_refine + i;  

                    if(0 == ind % FSLib->nx_refine) 
                    {
                        countX = 0;

                        if(0 != ind) 
                        {
                            countY++; 
                        }
                    }

                    ind2 = countY * FSLib->nx_fs + countX;

                    topo_fs[ind2] = topo_refine[ind]; 

                    countX++;                    
                }
            }

            ierr = VecRestoreArray(FSLib->gtopo_refine,  &topo_refine);    CHKERRQ(ierr);
        } 
    }
    else
    {
        ierr = VecGetArray(FSLib->gtopo_nug, &topo_nug);    CHKERRQ(ierr);

        if(1 == FSLib->refine)
        {
            // save value
            for(j = 0; j < fsY->nodes; j++)
            {
                for(i = 0; i < fsX->nodes; i++)
                {      
                    ind = j * fsX->nodes + i;

                    topo_nug[ind] =  topo_pass_f[ind] / scal->length_fs;

                    if(topo_nug[ind] > FSLib->rangeZ_end)   topo_nug[ind] = FSLib->rangeZ_end;
                    if(topo_nug[ind] < FSLib->rangeZ_begin) topo_nug[ind] = FSLib->rangeZ_begin;     
                }
            }

            // interploate
            for(j = 0; j < FSLib->ny_fs; j++)
            {
                for(i = 0; i < FSLib->nx_fs; i++)
                {     
                    x_coor = FSLib->ncoor_ori_x[i];
                    y_coor = FSLib->ncoor_ori_y[j];

                    dx     = fsX->dx;
                    dy     = fsY->dx;

                    // get nearest four index
                    // x-direction
                    m  = floor( (x_coor - x_begin) / dx );
                    mm = m + 1;
                    if( fsX->nodes == mm) mm -= 1;

                    // y-direction
                    n  = floor( (y_coor - y_begin) / dy );
                    nn = n + 1;
                    if( fsY->nodes == nn) nn -= 1;

                    // interpolate
                    ind   = j  * FSLib->nx_fs + i;    
                    ind_a = n  * fsX->nodes   + m;
                    ind_b = n  * fsX->nodes   + mm;
                    ind_c = nn * fsX->nodes   + m;
                    ind_d = nn * fsX->nodes   + mm;                             
                    // bilinear interpolation
                    wtx   = x_coor - fsX->ncoor[m];
                    wty   = y_coor - fsY->ncoor[n];
                                                                
                    // boundary (bottom or right)
                    if(m == mm)
                    {
                        if(n != nn)
                        {
                            topo_fs[ind] = ( (dy - wty) / dy )  * topo_nug[ind_a]+ \
                                                ( wty / dy )  * topo_nug[ind_c];  
                        }                     
                    }
                    else if(n == nn)
                    {
                        if(m != mm)
                        {
                            topo_fs[ind] = ( (dx - wtx) / dx ) * topo_nug[ind_a]+ \
                                                ( wtx / dx ) * topo_nug[ind_b];
                        }
                        else
                        {
                            topo_fs[ind] = topo_nug[ind_a];
                        }
                    }
                    else
                    {
                        topo_fs[ind] = ( (dy - wty) / dy ) * ( (dx - wtx) / dx ) * topo_nug[ind_a]+ \
                                            ( (dy - wty) / dy ) * ( wtx / dx ) * topo_nug[ind_b] + \
                                            ( wty / dy ) * ( (dx - wtx) / dx ) * topo_nug[ind_c] + \
                                            ( wty / dy ) * ( wtx / dx ) * topo_nug[ind_d]; 
                    }
                }
            }
        }
        else
        {
            PetscInt countX = 0, countY = 0;

            ierr = VecGetArray(FSLib->gtopo_refine, &topo_refine);         CHKERRQ(ierr);

            // save value
            for(j = 0; j < FSLib->ny_refine; j++)
            {
                for (i = 0; i < FSLib->nx_refine; i++)
                {
                    ind = j * FSLib->nx_refine + i;

                    topo_refine[ind] = topo_pass_f[ind] / scal->length_fs; // m(FastScape) to km(LaMEM)

                    if(topo_refine[ind] > FSLib->rangeZ_end)   topo_refine[ind] = FSLib->rangeZ_end;
                    if(topo_refine[ind] < FSLib->rangeZ_begin) topo_refine[ind] = FSLib->rangeZ_begin;
                }
            }

            for(j = 0; j < FSLib->ny_refine; j += FSLib->refine) 
            {   
                for(i = 0; i < FSLib->nx_refine; i += FSLib->refine) 
                {
                    ind = j * FSLib->nx_refine + i;  

                    if(0 == ind % FSLib->nx_refine) 
                    {
                        countX = 0;

                        if(0 != ind) 
                        {
                            countY++; 
                        }
                    }

                    ind2 = countY * fsX->nodes + countX;

                    topo_nug[ind2] = topo_refine[ind]; 

                    countX++;                    
                }
            }

            // interploate
            for(j = 0; j < FSLib->ny_fs; j++)
            {
                for(i = 0; i < FSLib->nx_fs; i++)
                {     
                    x_coor = FSLib->ncoor_ori_x[i];
                    y_coor = FSLib->ncoor_ori_y[j];

                    dx     = fsX->dx;
                    dy     = fsY->dx;

                    // get nearest four index
                    // x-direction
                    m  = floor( (x_coor - x_begin) / dx );
                    mm = m + 1;
                    if( fsX->nodes == mm) mm -= 1;

                    // y-direction
                    n  = floor( (y_coor - y_begin) / dy );
                    nn = n + 1;
                    if( fsY->nodes == nn) nn -= 1;

                    // interpolate
                    ind   = j  * FSLib->nx_fs + i;    
                    ind_a = n  * fsX->nodes   + m;
                    ind_b = n  * fsX->nodes   + mm;
                    ind_c = nn * fsX->nodes   + m;
                    ind_d = nn * fsX->nodes   + mm;                             
                    // bilinear interpolation
                    wtx   = x_coor - fsX->ncoor[m];
                    wty   = y_coor - fsY->ncoor[n];
                                                                
                    // boundary (bottom or right)
                    if(m == mm)
                    {
                        if(n != nn)
                        {
                            topo_fs[ind] = ( (dy - wty) / dy )  * topo_nug[ind_a]+ \
                                                ( wty / dy )  * topo_nug[ind_c];  
                        }                     
                    }
                    else if(n == nn)
                    {
                        if(m != mm)
                        {
                            topo_fs[ind] = ( (dx - wtx) / dx ) * topo_nug[ind_a]+ \
                                                ( wtx / dx ) * topo_nug[ind_b];
                        }
                        else
                        {
                            topo_fs[ind] = topo_nug[ind_a];
                        }
                    }
                    else
                    {
                        topo_fs[ind] = ( (dy - wty) / dy ) * ( (dx - wtx) / dx ) * topo_nug[ind_a]+ \
                                            ( (dy - wty) / dy ) * ( wtx / dx ) * topo_nug[ind_b] + \
                                            ( wty / dy ) * ( (dx - wtx) / dx ) * topo_nug[ind_c] + \
                                            ( wty / dy ) * ( wtx / dx ) * topo_nug[ind_d]; 
                    }                            
                }
            }

            ierr = VecRestoreArray(FSLib->gtopo_refine,  &topo_refine);    CHKERRQ(ierr);
        } 
    
        ierr = VecRestoreArray(FSLib->gtopo_nug, &topo_nug);    CHKERRQ(ierr);
    }

    PetscFunctionReturn(0); 
}
//---------------------------------------------------------------------------
PetscErrorCode PVSurfWriteVTSFS(FastScapeLib *FSLib, const char *dirName, PetscScalar *topo, PetscInt mode)
{
    FILE      *fp;
    Scaling   *scal;
    char      *fname;
    FreeSurf  *surf;
    size_t    offset = 0;
    PVSurf    *pvsurf;
    PetscErrorCode ierr;
    PetscFunctionBeginUser;

    PetscInt nx, ny;
    PetscScalar *silt_fraction = PETSC_NULL, *basement = PETSC_NULL, *total_erosion = PETSC_NULL;
    PetscScalar *drainage_area = PETSC_NULL, *erosion_rate = PETSC_NULL, *slope = PETSC_NULL;
    PetscScalar *curvature = PETSC_NULL, *chi = PETSC_NULL, *catchment = PETSC_NULL, *lake_depth = PETSC_NULL;

    // access context
    pvsurf = FSLib->pvsurf;
    surf = pvsurf->surf;
    scal = surf->jr->scal;
    
    // only processor 0 run the code
    if(!ISRankZero(PETSC_COMM_WORLD)) PetscFunctionReturn(0);

    fp = NULL;
    
    // initialize
    nx = FSLib->nx_solve;
    ny = FSLib->ny_solve;

    // open outfile_p_XXXXXX.vts file in the output directory (write mode)
    asprintf(&fname, "%s/%s_p0.vts", dirName, FSLib->outfile_fs);
 
    fp = fopen(fname,"wb");
    if(fp == NULL) SETERRQ(PETSC_COMM_SELF, 1,"cannot open file %s", fname);
    free(fname);

    // write header
    WriteXMLHeader(fp, "StructuredGrid");

    // open structured grid data block (write total grid size)
    fprintf(fp, "\t<StructuredGrid WholeExtent=\"%lld %lld %lld %lld 1 1\">\n",
            (LLD)(1), (LLD)(nx),
            (LLD)(1), (LLD)(ny));

    // open sub-domain (piece) description block
    fprintf(fp, "\t\t<Piece Extent=\"%lld %lld %lld %lld 1 1\">\n",
            (LLD)(1), (LLD)(nx),
            (LLD)(1), (LLD)(ny));

    // write cell data block (empty)
    fprintf(fp, "\t\t\t<CellData>\n");
    fprintf(fp, "\t\t\t</CellData>\n");

    // write coordinate block
    fprintf(fp, "\t\t<Points>\n");

    fprintf(fp,"\t\t\t<DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%lld\"/>\n",
            (LLD)offset);

    offset += sizeof(uint64_t) + sizeof(float)*(size_t)(nx * ny * 3);

    fprintf(fp, "\t\t</Points>\n");

    // write description of output vectors
    fprintf(fp, "\t\t<PointData>\n");

    // fastscape grid
    // output topo
    if(FSLib->out_topofs)
    {
        fprintf(fp,"\t\t\t<DataArray type=\"Float32\" Name=\"topoFs %s\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n",
            scal->lbl_length, (LLD)offset);

        offset += sizeof(uint64_t) + sizeof(float)*(size_t)(nx * ny);
    }
    
    // output silt fraction
    if(FSLib->out_silt_fraction)
    {
        fprintf(fp,"\t\t\t<DataArray type=\"Float32\" Name=\"silt_fraction %s\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n",
            scal->lbl_unit, (LLD)offset);

        offset += sizeof(uint64_t) + sizeof(float)*(size_t)(nx * ny);
    }
    // output basement
    if(FSLib->out_basement)
    {
        fprintf(fp,"\t\t\t<DataArray type=\"Float32\" Name=\"basement %s\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n",
            scal->lbl_length, (LLD)offset);

        offset += sizeof(uint64_t) + sizeof(float)*(size_t)(nx * ny);
    }
    // output total_erosion
    if(FSLib->out_total_erosion)
    {
        fprintf(fp,"\t\t\t<DataArray type=\"Float32\" Name=\"total_erosion %s\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n",
           scal->lbl_length, (LLD)offset);

        offset += sizeof(uint64_t) + sizeof(float)*(size_t)(nx * ny);
    }
    // output drainage_area
    if(FSLib->out_drainage_area)
    {
        fprintf(fp,"\t\t\t<DataArray type=\"Float32\" Name=\"drainage_area %s\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n",
            scal->lbl_area_fs, (LLD)offset);

        offset += sizeof(uint64_t) + sizeof(float)*(size_t)(nx * ny);
    }
    // output erosion_rate
    if(FSLib->out_erosion_rate)
    {
        fprintf(fp,"\t\t\t<DataArray type=\"Float32\" Name=\"erosion_rate %s\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n",
            scal->lbl_rate, (LLD)offset);

        offset += sizeof(uint64_t) + sizeof(float)*(size_t)(nx * ny);
    }
    // output slope
    if(FSLib->out_slope)
    {
        fprintf(fp,"\t\t\t<DataArray type=\"Float32\" Name=\"slope %s\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n",
            scal->lbl_degree, (LLD)offset);

        offset += sizeof(uint64_t) + sizeof(float)*(size_t)(nx * ny);
    }
    // output curvature
    if(FSLib->out_curvature)
    {
        fprintf(fp,"\t\t\t<DataArray type=\"Float32\" Name=\"curvature %s\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n",
            scal->lbl_unit, (LLD)offset);

        offset += sizeof(uint64_t) + sizeof(float)*(size_t)(nx * ny);
    }
    // output chi
    if(FSLib->out_chi)
    {
        fprintf(fp,"\t\t\t<DataArray type=\"Float32\" Name=\"chi %s\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n",
            scal->lbl_length, (LLD)offset);

        offset += sizeof(uint64_t) + sizeof(float)*(size_t)(nx * ny);
    }
    // output catchment
    if(FSLib->out_catchment)
    {
        fprintf(fp,"\t\t\t<DataArray type=\"Float32\" Name=\"catchment %s\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n",
            scal->lbl_area_fs, (LLD)offset);

        offset += sizeof(uint64_t) + sizeof(float)*(size_t)(nx * ny);
    }
    // output lake_depth
    if(FSLib->out_lake_depth)
    {
        fprintf(fp,"\t\t\t<DataArray type=\"Float32\" Name=\"lake_depth %s\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n",
            scal->lbl_length, (LLD)offset);

        offset += sizeof(uint64_t) + sizeof(float)*(size_t)(nx * ny);
    }   
    
    fprintf(fp, "\t\t</PointData>\n");

    // close sub-domain and grid blocks
    fprintf(fp, "\t\t</Piece>\n");
    fprintf(fp, "\t</StructuredGrid>\n");

    // write appended data section
    fprintf(fp, "\t<AppendedData encoding=\"raw\">\n");
    fprintf(fp,"_");
    
    // write point coordinates
    // allocate output buffer
    ierr = PetscMalloc((size_t)(_max_num_comp_surf_ * nx * ny)*sizeof(float), &FSLib->buff_fs);    CHKERRQ(ierr);

    ierr = PVSurfWriteCoordFS (FSLib, fp, topo, mode); CHKERRQ(ierr);

// topography
    ierr = PVSurfWriteInfFS  (FSLib, fp, topo, 1);           CHKERRQ(ierr);
    // silt fraction
    if(FSLib->out_silt_fraction) 
    {
        ierr = VecGetArray(FSLib->silt_fraction, &silt_fraction);  CHKERRQ(ierr);
        ierr = PVSurfWriteInfFS  (FSLib, fp, silt_fraction, 2);       CHKERRQ(ierr);
        ierr = VecRestoreArray(FSLib->silt_fraction, &silt_fraction);  CHKERRQ(ierr);
    }
    // basement
    if(FSLib->out_basement) 
    {
        ierr = VecGetArray(FSLib->basement, &basement);  CHKERRQ(ierr);
        ierr = PVSurfWriteInfFS  (FSLib, fp, basement, 3);       CHKERRQ(ierr);
        ierr = VecRestoreArray(FSLib->basement, &basement);  CHKERRQ(ierr);
    }
    // total_erosion
    if(FSLib->out_total_erosion) 
    {
        ierr = VecGetArray(FSLib->total_erosion, &total_erosion);  CHKERRQ(ierr);
        ierr = PVSurfWriteInfFS  (FSLib, fp, total_erosion, 4);       CHKERRQ(ierr);
        ierr = VecRestoreArray(FSLib->total_erosion, &total_erosion);  CHKERRQ(ierr);
    }
    // drainage_area
    if(FSLib->out_drainage_area) 
    {
        ierr = VecGetArray(FSLib->drainage_area, &drainage_area);  CHKERRQ(ierr);
        ierr = PVSurfWriteInfFS  (FSLib, fp, drainage_area, 5);       CHKERRQ(ierr);
        ierr = VecRestoreArray(FSLib->drainage_area, &drainage_area);  CHKERRQ(ierr);
    }
    // erosion_rate
    if(FSLib->out_erosion_rate) 
    {
        ierr = VecGetArray(FSLib->erosion_rate, &erosion_rate);  CHKERRQ(ierr);
        ierr = PVSurfWriteInfFS  (FSLib, fp, erosion_rate, 6);       CHKERRQ(ierr);
        ierr = VecRestoreArray(FSLib->erosion_rate, &erosion_rate);  CHKERRQ(ierr);
    }
    // slope
    if(FSLib->out_slope) 
    {
        ierr = VecGetArray(FSLib->slope, &slope);  CHKERRQ(ierr);
        ierr = PVSurfWriteInfFS  (FSLib, fp, slope, 7);       CHKERRQ(ierr);
        ierr = VecRestoreArray(FSLib->slope, &slope);  CHKERRQ(ierr);
    }
    // curvature
    if(FSLib->out_curvature) 
    {
        ierr = VecGetArray(FSLib->curvature, &curvature);  CHKERRQ(ierr);
        ierr = PVSurfWriteInfFS  (FSLib, fp, curvature, 8);       CHKERRQ(ierr);
        ierr = VecRestoreArray(FSLib->curvature, &curvature);  CHKERRQ(ierr);
    }
    // chi
    if(FSLib->out_chi) 
    {
        ierr = VecGetArray(FSLib->chi, &chi);  CHKERRQ(ierr);
        ierr = PVSurfWriteInfFS  (FSLib, fp, chi, 9);       CHKERRQ(ierr);
        ierr = VecRestoreArray(FSLib->chi, &chi);  CHKERRQ(ierr);
    }
    // catchment
    if(FSLib->out_catchment) 
    {
        ierr = VecGetArray(FSLib->catchment, &catchment);  CHKERRQ(ierr);
        ierr = PVSurfWriteInfFS  (FSLib, fp, catchment, 10);       CHKERRQ(ierr);
        ierr = VecRestoreArray(FSLib->catchment, &catchment);  CHKERRQ(ierr);
    }
    // lake_depth
    if(FSLib->out_lake_depth) 
    {
        ierr = VecGetArray(FSLib->lake_depth, &lake_depth);  CHKERRQ(ierr);
        ierr = PVSurfWriteInfFS  (FSLib, fp, lake_depth, 11);       CHKERRQ(ierr);
        ierr = VecRestoreArray(FSLib->lake_depth, &lake_depth);  CHKERRQ(ierr);
    }
    // close appended data section and file
    fprintf(fp, "\n\t</AppendedData>\n");
    fprintf(fp, "</VTKFile>\n");

    // close file
    fclose(fp);

    ierr = PetscFree(FSLib->buff_fs);    CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PVSurfWriteCoordFS(FastScapeLib *FSLib, FILE *fp, PetscScalar *topo, PetscInt mode)
{
    float       *buff;
    PetscInt    i, j, ind, cn, nx, ny;
    FSGrid  *fsX;
    FSGrid  *fsY;
    PetscFunctionBeginUser;

    if(!ISRankZero(PETSC_COMM_WORLD)) PetscFunctionReturn(0);

    fsX = &FSLib->fsX;
    fsY = &FSLib->fsY;

    // initialize
    buff   = FSLib->buff_fs;
    
    nx = FSLib->nx_solve;
    ny = FSLib->ny_solve;

    cn     = 0;

    for(j = 0; j < ny; j++)
    {
        for(i = 0; i < nx; i++)
        {
            ind = j * nx + i;

            if( 0 == mode )
            {
                // store node coordinates
                buff[cn++] = (float)(fsX->ncoor[i]); 
                buff[cn++] = (float)(fsY->ncoor[j]);
                buff[cn++] = (float)(topo[ind] * FSLib->vec_times); // km -> m
            }
            else if( 2 == mode )
            {
                // store node coordinates
                buff[cn++] = (float)(fsX->ncoor_extend[i]); 
                buff[cn++] = (float)(fsY->ncoor_extend[j]);
                buff[cn++] = (float)(topo[ind] * FSLib->vec_times); // km -> m                
            }
            else if(1 == mode || 3 == mode)
            {
                buff[cn++] = (float)(fsX->ncoor_refine[i]); 
                buff[cn++] = (float)(fsY->ncoor_refine[j]);
                buff[cn++] = (float)(topo[ind] * FSLib->vec_times); // km -> m                
            }
        }
    }

    OutputBufferWrite(fp, buff, cn);

    PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PVSurfWriteInfFS(FastScapeLib *FSLib, FILE *fp, PetscScalar *Inf, PetscInt InfMode)
{
    float       *buff;
    PetscInt    i, j, ind, cn, nx, ny;
    Scaling *scal;
    PetscFunctionBeginUser;

    scal = FSLib->scal;

    // initialize
    buff   = FSLib->buff_fs;

    nx = FSLib->nx_solve;
    ny = FSLib->ny_solve;

    cn   = 0;

    for(j = 0; j < ny; j++)
    {
        for(i = 0; i < nx; i++)
        {
            ind = j * nx + i;
            // topography
            if(1 == InfMode)
            {
                buff[cn++] = (float)(Inf[ind]);
            }
            // silt fraction
            if(2 == InfMode)
            {
                buff[cn++] = (float)(Inf[ind]);
            }
            // basement
            if(3 == InfMode)
            {
                buff[cn++] = (float)(Inf[ind]/scal->length_fs);
            }
            // total_erosion
            if(4 == InfMode)
            {
                buff[cn++] = (float)(Inf[ind]/scal->length_fs);
            }
            // drainage_area
            if(5 == InfMode)
            {
                buff[cn++] = (float)(Inf[ind]/scal->area_fs);
            }
            // erosion
            if(6 == InfMode)
            {
                buff[cn++] = (float)(Inf[ind]/scal->length_fs);
            }
            // slope
            if(7 == InfMode)
            {
                buff[cn++] = (float)(Inf[ind]);
            }
            // curvature
            if(8 == InfMode)
            {
                buff[cn++] = (float)(Inf[ind]);
            }
            // chi
            if(9 == InfMode)
            {
                buff[cn++] = (float)(Inf[ind]/scal->length_fs);
            }
            // catchment
            if(10 == InfMode)
            {
                buff[cn++] = (float)(Inf[ind]/scal->area_fs);
            }
            // lake_depth
            if(11 == InfMode)
            {
                buff[cn++] = (float)(Inf[ind]/scal->length_fs);
            }
        }
    }

    OutputBufferWrite(fp, buff, cn);
    
    PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode UpdatePVDFileFS(
        const char *dirName, const char *outfile, const char *ext,
        long int *offset, PetscScalar ttime, PetscInt outpvd, PetscInt step)
{
    FILE        *fp;
    char        *fname;

    PetscErrorCode ierr;
    PetscFunctionBeginUser;

    // check whether pvd is requested
    if(!outpvd) PetscFunctionReturn(0);

    // only first process generates this file (WARNING! Bottleneck!)
    if(!ISRankZero(PETSC_COMM_WORLD)) PetscFunctionReturn(0);

    // open outfile.pvd file (write or update mode)
    asprintf(&fname, "%s.pvd", outfile);
    if(step == 1) fp = fopen(fname,"wb");
    else       fp = fopen(fname,"r+b");
    free(fname);

    if(fp == NULL) SETERRQ(PETSC_COMM_SELF, 1,"cannot open file %s", fname);

    if(step == 1)
    {
        // write header
        WriteXMLHeader(fp, "Collection");

        // open time step collection
        fprintf(fp,"<Collection>\n");
    }
    else
    {
        // put the file pointer on the next entry
        ierr = fseek(fp, (*offset), SEEK_SET); CHKERRQ(ierr);
    }

    // add entry to .pvd file
    fprintf(fp,"\t<DataSet timestep=\"%1.6e\" file=\"%s/%s_%s\"/>\n",
        ttime, dirName, outfile, ext);

    // store current position in the file
    (*offset) = ftell(fp);

    // close time step collection
    fprintf(fp,"</Collection>\n");
    fprintf(fp,"</VTKFile>\n");

    // close file
    fclose(fp);

    PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode SavePvtsFS(FastScapeLib *FSLib, PetscScalar ttime, PetscInt step, const char *dirName, PetscScalar *topo)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    PetscInt mode;

    // check activation
    if(!FSLib->outsurf_fs) PetscFunctionReturn(0);

    // update .pvd file if necessary
    ierr = UpdatePVDFileFS(dirName, FSLib->outfile_fs, "p0.vts", &FSLib->offset_fs, ttime, FSLib->outpvd_fs, step); CHKERRQ(ierr);

    // write sub-domain data .vts files
    if(0 == FSLib->fs2D)
    {
        if ( 1 == FSLib->refine) mode = 0;
        else mode = 1;
    }
    else
    {
        if ( 1 == FSLib->refine) mode = 2;
        else mode = 3;        
    }

    ierr = PVSurfWriteVTSFS(FSLib, dirName, topo, mode);           CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FastScapeSave(FastScapeLib *FSLib, PetscInt step_fs, PetscScalar time_fs)
{
    PetscErrorCode ierr;
    PetscInt status;
    PetscScalar *gtopo_fs = PETSC_NULL, *gtopo_extend = PETSC_NULL, *topo_refine = PETSC_NULL, *topo_et_refine = PETSC_NULL;
    char        *dirName;
    
    // create directory(encode current time & steo number)    
    // update time stamp and counter
    step_fs++;

    asprintf(&dirName, "Timestep_%1.8lld_%1.8e", (LLD)step_fs, time_fs);

    // create output directory
    #ifdef _WIN32
    // call this on windows machines
    status = mkdir(dirName);
    #else
    // standard access pattern drwxr-xr-x
    status = mkdir(dirName, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
    #endif 
    if(status && errno != EEXIST)
    {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Failed to create directory %s", dirName);
    }
    
    // only saved in processor 0
    if(ISRankZero(PETSC_COMM_WORLD))
    {
        
        if(FSLib->outsurf_fs)
        {
            if( 0 == FSLib->fs2D)
            {
                if(1 == FSLib->refine)
                {
                    if( 0 == FSLib->non_uniform_grid )
                    {
                        ierr = VecGetArray(FSLib->gtopo_fs,  &gtopo_fs);                                 CHKERRQ(ierr);
                        ierr = SavePvtsFS(FSLib, time_fs, step_fs, dirName, gtopo_fs);                   CHKERRQ(ierr);  
                        ierr = VecRestoreArray(FSLib->gtopo_fs,  &gtopo_fs);                             CHKERRQ(ierr);      
                    }
                    else
                    {
                        ierr = VecGetArray(FSLib->gtopo_nug,  &gtopo_fs);                                 CHKERRQ(ierr);
                        ierr = SavePvtsFS(FSLib, time_fs, step_fs, dirName, gtopo_fs);                   CHKERRQ(ierr);  
                        ierr = VecRestoreArray(FSLib->gtopo_nug,  &gtopo_fs);                             CHKERRQ(ierr); 
                    }
                }
                else
                {
                    ierr = VecGetArray(FSLib->gtopo_refine,  &topo_refine);                          CHKERRQ(ierr);
                    ierr = SavePvtsFS(FSLib, time_fs, step_fs, dirName, topo_refine);          CHKERRQ(ierr);
                    ierr = VecRestoreArray(FSLib->gtopo_refine,  &topo_refine);                      CHKERRQ(ierr);
                }
            }
            else
            {
                if(1 == FSLib->refine)
                {
                    ierr = VecGetArray(FSLib->gtopo_extend,  &gtopo_extend);                         CHKERRQ(ierr);
                    ierr = SavePvtsFS(FSLib, time_fs, step_fs, dirName, gtopo_extend);         CHKERRQ(ierr);  
                    ierr = VecRestoreArray(FSLib->gtopo_extend,  &gtopo_extend);                     CHKERRQ(ierr); 
                }
                else
                {
                    ierr = VecGetArray(FSLib->gtopo_et_refine,  &topo_et_refine);                     CHKERRQ(ierr);
                    ierr = SavePvtsFS(FSLib, time_fs, step_fs, dirName, topo_et_refine);  CHKERRQ(ierr);  
                    ierr = VecRestoreArray(FSLib->gtopo_et_refine,  &topo_et_refine);                 CHKERRQ(ierr);                     
                }
            }  
        }   
    }

    PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FastScapeFortranCppAdvc(FastScapeLib *FSLib, PetscScalar dt_max, PetscScalar dt_n, PetscInt nstep, 
    PetscScalar step_fs, PetscScalar *vx_pass, PetscScalar *vy_pass, PetscScalar *vz_pass, PetscScalar *topo_pass)
{
    PetscErrorCode ierr;
    PetscInt istep, ind, i, j;
    PetscScalar *topo_random = PETSC_NULL, *kf = PETSC_NULL, *kd = PETSC_NULL;
    PetscScalar *silt_fraction = PETSC_NULL, *basement = PETSC_NULL, *total_erosion = PETSC_NULL;
    PetscScalar *drainage_area = PETSC_NULL, *erosion_rate = PETSC_NULL, *slope = PETSC_NULL;
    PetscScalar *curvature = PETSC_NULL, *chi = PETSC_NULL, *catchment = PETSC_NULL, *lake_depth = PETSC_NULL;
    char *endptr; 
    PetscInt ibc_int;

    // initialize FastScape
    fastscape_init_();

    fastscape_set_nx_ny_(&FSLib->nx_solve, &FSLib->ny_solve);    

    if(1 == FSLib->random_noise)
    {
        topo_random = (PetscScalar *)malloc(FSLib->nodes_solve * sizeof(PetscScalar));        
    }

    kf          = (PetscScalar *)malloc(FSLib->nodes_solve * sizeof(PetscScalar)); 
    kd          = (PetscScalar *)malloc(FSLib->nodes_solve * sizeof(PetscScalar)); 

    // output
    if(FSLib->out_silt_fraction)
    {
        ierr = VecGetArray(FSLib->silt_fraction, &silt_fraction);  CHKERRQ(ierr);
    }
    // output basement
    if(FSLib->out_basement)
    {
        ierr = VecGetArray(FSLib->basement,      &basement);       CHKERRQ(ierr);
    }
    // output total_erosion
    if(FSLib->out_total_erosion)
    {
        ierr = VecGetArray(FSLib->total_erosion, &total_erosion);  CHKERRQ(ierr);
    }
    // output drainage_area
    if(FSLib->out_drainage_area)
    {
        ierr = VecGetArray(FSLib->drainage_area, &drainage_area);  CHKERRQ(ierr);
    }
    // output erosion_rate
    if(FSLib->out_erosion_rate)
    {
        ierr = VecGetArray(FSLib->erosion_rate, &erosion_rate);   CHKERRQ(ierr);
    }
    // output slope
    if(FSLib->out_slope)
    {
        ierr = VecGetArray(FSLib->slope,       &slope);           CHKERRQ(ierr);
    }
    // output curvature
    if(FSLib->out_curvature)
    {
        ierr = VecGetArray(FSLib->curvature,   &curvature);       CHKERRQ(ierr);
    }
    // output chi
    if(FSLib->out_chi)
    {
        ierr = VecGetArray(FSLib->chi,       &chi);               CHKERRQ(ierr);
    }
    // output catchment
    if(FSLib->out_catchment)
    {
        ierr = VecGetArray(FSLib->catchment, &catchment);         CHKERRQ(ierr);
    }
    // output lake_depth
    if(FSLib->out_lake_depth)
    {
        ierr = VecGetArray(FSLib->lake_depth, &lake_depth);       CHKERRQ(ierr);
    }     

    // allocate memory
    fastscape_setup_();

    // set model dimensions
    if(0 == FSLib->fs2D)
    {
        fastscape_set_xl_yl_(&FSLib->rangeX, &FSLib->rangeY);
    }
    else
    {
        fastscape_set_xl_yl_(&FSLib->extendedXRange, &FSLib->extendedYRange);
    }

    // set time step
    fastscape_set_dt_(&dt_max);

    // random noise (uniform distribution)
    mt19937 generator;
    uniform_int_distribution<int> distribution(1, 10000);

    // set random initial topography & erosional parameters & vz boundary condition
    if(0 == step_fs) 
    {
        for(j = 0; j < FSLib->ny_solve; j++)
        {
            for(i = 0; i < FSLib->nx_solve; i++)
            {
                ind = j * FSLib->nx_solve + i;

                if(1 == FSLib->random_noise)
                {
                    topo_random[ind] =  distribution(generator) / 10000.0;
                    topo_pass[ind]  += topo_random[ind];
                }

                kf[ind]          = FSLib->kf;
                kd[ind]          = FSLib->kd;


                // i for nx; j for ny
                // bottom (j==0) -- right(i==end) -- top (j==end) -- left (i==0)
                if( '1' == FSLib->FS_VELBC[0] ) 
                {
                    if(  0 == j )
                    {
                        vx_pass[ind] = 0;
                        vy_pass[ind] = 0;
                        vz_pass[ind] = 0;                  
                    }
                }
                if( '1' == FSLib->FS_VELBC[1] ) 
                {
                    if( (FSLib->nx_solve - 1) == i ) 
                    {
                        vx_pass[ind] = 0;
                        vy_pass[ind] = 0;
                        vz_pass[ind] = 0;                  
                    }              
                }
                if( '1' == FSLib->FS_VELBC[2] ) 
                {
                    if( (FSLib->ny_solve - 1) == j )
                    {
                        vx_pass[ind] = 0;
                        vy_pass[ind] = 0;
                        vz_pass[ind] = 0;                  
                    }
                }
                if( '1' == FSLib->FS_VELBC[3] ) 
                {
                    if( 0 == i )
                    {
                        vx_pass[ind] = 0;
                        vy_pass[ind] = 0;
                        vz_pass[ind] = 0;                  
                    }
                }

            }
        } 
    } 
    else
    {
        for(j = 0; j < FSLib->ny_solve; j++)
        {
            for(i = 0; i < FSLib->nx_solve; i++)
            {
                ind = j * FSLib->nx_solve + i;

                kf[ind] = FSLib->kf;
                kd[ind] = FSLib->kd;

                // i for nx; j for ny
                // bottom (j==0) -- right(i==end) -- top (j==end) -- left (i==0)
                if( '1' == FSLib->FS_VELBC[0] ) 
                {
                    if(  0 == j )
                    {
                        vx_pass[ind] = 0;
                        vy_pass[ind] = 0;
                        vz_pass[ind] = 0;                  
                    }
                }
                if( '1' == FSLib->FS_VELBC[1] ) 
                {
                    if( (FSLib->nx_solve - 1) == i ) 
                    {
                        vx_pass[ind] = 0;
                        vy_pass[ind] = 0;
                        vz_pass[ind] = 0;                  
                    }              
                }
                if( '1' == FSLib->FS_VELBC[2] ) 
                {
                    if( (FSLib->ny_solve - 1) == j )
                    {
                        vx_pass[ind] = 0;
                        vy_pass[ind] = 0;
                        vz_pass[ind] = 0;                  
                    }
                }
                if( '1' == FSLib->FS_VELBC[3] ) 
                {
                    if( 0 == i )
                    {
                        vx_pass[ind] = 0;
                        vy_pass[ind] = 0;
                        vz_pass[ind] = 0;                  
                    }
                }

            }
        }         
    }

    fastscape_init_h_(topo_pass);

    // set erosional parameters
    fastscape_set_erosional_parameters_(kf, &FSLib->kfsed, &FSLib->m, &FSLib->n, kd, &FSLib->kdsed, &FSLib->g, &FSLib->gsed, &FSLib->p);

    // set marine transport parameters
    if (1 == FSLib->setMarine)
    {
        fastscape_set_marine_parameters_(&FSLib->sealevel, &FSLib->poro_silt, &FSLib->poro_sand, &FSLib->zporo_silt, 
                                         &FSLib->zporo_sand, &FSLib->ratio, &FSLib->Lsolve, &FSLib->kds_silt, &FSLib->kds_sand);
    }
    
    // set vz
    fastscape_set_u_(vz_pass);
   
    // set vx, vy
    fastscape_set_v_(vx_pass, vy_pass);

    ibc_int = strtol(FSLib->FS_BC, &endptr, 10);

    // set boundary conditions
    fastscape_set_bc_(&ibc_int);

    // set number of time steps and initialize counter istep
    fastscape_get_step_(&istep);

    // loop on time stepping
    do 
    {
        if( (nstep - 1) > istep)
        {
            // execute step
            fastscape_execute_step_();
            // get value of time step counter
            fastscape_get_step_(&istep);
        }
        else
        {
            if(0 < dt_n)
            {
                // reset dt
                fastscape_set_dt_(&dt_n);            
                // execute step
                fastscape_execute_step_();
                // get value of time step counter
                fastscape_get_step_(&istep);
            }
            else
            {
                // execute step
                fastscape_execute_step_();
                // get value of time step counter
                fastscape_get_step_(&istep);
            }
           // output
            // output h values
            fastscape_copy_h_(topo_pass);
            // output silt fraction
            if(FSLib->out_silt_fraction)
            {
                fastscape_copy_f_(silt_fraction);
            }
            // output basement
            if(FSLib->out_basement)
            {
                fastscape_copy_basement_(basement);
            }
            // output total_erosion
            if(FSLib->out_total_erosion)
            {
                fastscape_copy_total_erosion_(total_erosion);
            }
            // output drainage_area
            if(FSLib->out_drainage_area)
            {
                fastscape_copy_drainage_area_(drainage_area);
            }
            // output erosion_rate
            if(FSLib->out_erosion_rate)
            {
                fastscape_copy_erosion_rate_(erosion_rate);
            }
            // output slope
            if(FSLib->out_slope)
            {
                fastscape_copy_slope_(slope);
            }
            // output curvature
            if(FSLib->out_curvature)
            {
                fastscape_copy_curvature_(curvature);
            }
            // output chi
            if(FSLib->out_chi)
            {
                fastscape_copy_chi_(chi);
            }
            // output catchment
            if(FSLib->out_catchment)
            {
                fastscape_copy_catchment_(catchment);
            }
            // output lake_depth
            if(FSLib->out_lake_depth)
            {
                fastscape_copy_lake_depth_(lake_depth);
            }          
        }
    } while (istep < nstep);

    // output timing
    fastscape_debug_();

    // end FastScape run
    fastscape_destroy_();

    if(1 == FSLib->random_noise)
    {
        free(topo_random);
        topo_random = PETSC_NULL;
    }

    free(kf);
    free(kd);

    kf          = PETSC_NULL;
    kd          = PETSC_NULL;


    // output
    if(FSLib->out_silt_fraction)
    {
        ierr = VecRestoreArray(FSLib->silt_fraction, &silt_fraction);  CHKERRQ(ierr);
    }
    // output basement
    if(FSLib->out_basement)
    {
        ierr = VecRestoreArray(FSLib->basement,      &basement);       CHKERRQ(ierr);
    }
    // output total_erosion
    if(FSLib->out_total_erosion)
    {
        ierr = VecRestoreArray(FSLib->total_erosion, &total_erosion);  CHKERRQ(ierr);
    }
    // output drainage_area
    if(FSLib->out_drainage_area)
    {
        ierr = VecRestoreArray(FSLib->drainage_area, &drainage_area);  CHKERRQ(ierr);
    }
    // output erosion_rate
    if(FSLib->out_erosion_rate)
    {
        ierr = VecRestoreArray(FSLib->drainage_area, &erosion_rate);   CHKERRQ(ierr);
    }
    // output slope
    if(FSLib->out_slope)
    {
        ierr = VecRestoreArray(FSLib->slope,        &slope);          CHKERRQ(ierr);
    }
    // output curvature
    if(FSLib->out_curvature)
    {
        ierr = VecRestoreArray(FSLib->curvature,    &curvature);      CHKERRQ(ierr);
    }
    // output chi
    if(FSLib->out_chi)
    {
        ierr = VecRestoreArray(FSLib->chi,          &chi);            CHKERRQ(ierr);
    }
    // output catchment
    if(FSLib->out_catchment)
    {
        ierr = VecRestoreArray(FSLib->catchment,    &catchment);      CHKERRQ(ierr);
    }
    // output lake_depth
    if(FSLib->out_lake_depth)
    {
        ierr = VecRestoreArray(FSLib->lake_depth,   &lake_depth);     CHKERRQ(ierr);
    }   

    PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FastScapeReadRestart(FastScapeLib *FSLib, FILE *fp)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;

    if(ISRankZero(PETSC_COMM_WORLD))
    {
        ierr =  FastScapeCreateSurfaceGrid(FSLib, 2); CHKERRQ(ierr); 

        if( 0 == FSLib->fs2D && 1 == FSLib->refine  )
        {
            // read topography vector
            ierr = VecReadRestart(FSLib->gtopo_fs, fp); CHKERRQ(ierr);
        }
        if( 1 == FSLib->non_uniform_grid)
        {
            // read topography vector
            ierr = VecReadRestart(FSLib->gtopo_nug, fp); CHKERRQ(ierr);            
        }
        if( 0 == FSLib->fs2D && FSLib->refine > 1  )
        {
            // read topography vector
            ierr = VecReadRestart(FSLib->gtopo_refine, fp); CHKERRQ(ierr);
        }
        else if( 1 == FSLib->fs2D && 1 == FSLib->refine )
        {
            // read topography vector
            ierr = VecReadRestart(FSLib->gtopo_extend, fp); CHKERRQ(ierr);        
        }
        else if( 1 == FSLib->fs2D && FSLib->refine > 1)
        {
            // read topography vector
            ierr = VecReadRestart(FSLib->gtopo_et_refine, fp); CHKERRQ(ierr);          
        }
    }   

    PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FastScapeWriteRestart(FastScapeLib *FSLib, FILE *fp)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;

    if(ISRankZero(PETSC_COMM_WORLD))
    {

        if( 0 == FSLib->fs2D && 1 == FSLib->refine  )
        {
            // store topography vector
            ierr = VecWriteRestart(FSLib->gtopo_fs, fp); CHKERRQ(ierr);
        }
        if(1 == FSLib->non_uniform_grid)
        {
            // store topography vector
            ierr = VecWriteRestart(FSLib->gtopo_nug, fp); CHKERRQ(ierr);            
        }
        if( 0 == FSLib->fs2D && FSLib->refine > 1)
        {
            // store topography vector
            ierr = VecWriteRestart(FSLib->gtopo_refine, fp); CHKERRQ(ierr);
        }
        else if( (1 == FSLib->fs2D && 1 == FSLib->refine)  )
        {
            // store topography vector
            ierr = VecWriteRestart(FSLib->gtopo_extend, fp); CHKERRQ(ierr);        
        }
        else if( 1 == FSLib->fs2D && FSLib->refine > 1)
        {
            // store topography vector 
            ierr = VecWriteRestart(FSLib->gtopo_et_refine, fp); CHKERRQ(ierr);         
        }
    }

    PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode FastScapeDestroy(FastScapeLib *FSLib)
{
    PetscFunctionBeginUser;
    PetscErrorCode ierr;

    ierr = VecDestroy(&FSLib->gtopo_nug); CHKERRQ(ierr);         
    ierr = VecDestroy(&FSLib->vx_nug); CHKERRQ(ierr);   
    ierr = VecDestroy(&FSLib->vy_nug); CHKERRQ(ierr);   
    ierr = VecDestroy(&FSLib->vz_nug); CHKERRQ(ierr);   
    ierr = VecDestroy(&FSLib->vx_fs); CHKERRQ(ierr);   
    ierr = VecDestroy(&FSLib->vy_fs); CHKERRQ(ierr);   
    ierr = VecDestroy(&FSLib->vz_fs); CHKERRQ(ierr);   
    ierr = VecDestroy(&FSLib->gtopo_fs); CHKERRQ(ierr);   
    ierr = VecDestroy(&FSLib->vx_collect); CHKERRQ(ierr);   
    ierr = VecDestroy(&FSLib->vy_collect); CHKERRQ(ierr);   
    ierr = VecDestroy(&FSLib->vz_collect); CHKERRQ(ierr);   
    ierr = VecDestroy(&FSLib->gtopo_collect); CHKERRQ(ierr);      
    ierr = VecDestroy(&FSLib->gtopo_refine); CHKERRQ(ierr);   
    ierr = VecDestroy(&FSLib->vx_refine); CHKERRQ(ierr);   
    ierr = VecDestroy(&FSLib->vy_refine); CHKERRQ(ierr);  
    ierr = VecDestroy(&FSLib->vz_refine); CHKERRQ(ierr);  
    ierr = VecDestroy(&FSLib->gtopo_extend); CHKERRQ(ierr);  
    ierr = VecDestroy(&FSLib->vx_extend); CHKERRQ(ierr); 
    ierr = VecDestroy(&FSLib->vy_extend); CHKERRQ(ierr); 
    ierr = VecDestroy(&FSLib->vz_extend); CHKERRQ(ierr); 
    ierr = VecDestroy(&FSLib->gtopo_et_refine); CHKERRQ(ierr); 
    ierr = VecDestroy(&FSLib->vx_et_refine); CHKERRQ(ierr);
    ierr = VecDestroy(&FSLib->vy_et_refine); CHKERRQ(ierr);
    ierr = VecDestroy(&FSLib->vz_et_refine); CHKERRQ(ierr);
    ierr = VecDestroy(&FSLib->silt_fraction); CHKERRQ(ierr);
    ierr = VecDestroy(&FSLib->basement); CHKERRQ(ierr);
    ierr = VecDestroy(&FSLib->total_erosion); CHKERRQ(ierr);
    ierr = VecDestroy(&FSLib->drainage_area); CHKERRQ(ierr);
    ierr = VecDestroy(&FSLib->erosion_rate); CHKERRQ(ierr);
    ierr = VecDestroy(&FSLib->slope); CHKERRQ(ierr);
    ierr = VecDestroy(&FSLib->curvature); CHKERRQ(ierr);
    ierr = VecDestroy(&FSLib->chi); CHKERRQ(ierr);
    ierr = VecDestroy(&FSLib->catchment); CHKERRQ(ierr);
    ierr = VecDestroy(&FSLib->lake_depth); CHKERRQ(ierr);

    PetscFree(FSLib->buff_fs);

    if(ISRankZero(PETSC_COMM_WORLD))
    {
        free(FSLib->fsX.ncoor);
        free(FSLib->fsY.ncoor);
        FSLib->fsX.ncoor = PETSC_NULL;
        FSLib->fsY.ncoor = PETSC_NULL;

        if(1 == FSLib->non_uniform_grid)
        {
            if(1 == FSLib->extendedX)
            {
                free(FSLib->ncoor_ori_y);
                FSLib->ncoor_ori_y = PETSC_NULL;
            }
            else
            {
                free(FSLib->ncoor_ori_x);
                FSLib->ncoor_ori_x = PETSC_NULL;
            }       
        }

        if(1 == FSLib->fs2D)
        {
            free(FSLib->fsX.ncoor_extend);
            free(FSLib->fsY.ncoor_extend);     
            FSLib->fsX.ncoor_extend = PETSC_NULL;
            FSLib->fsY.ncoor_extend = PETSC_NULL;
        }

        if(1 < FSLib->refine)
        {
            free(FSLib->fsX.ncoor_refine);
            free(FSLib->fsY.ncoor_refine);
            FSLib->fsX.ncoor_refine = PETSC_NULL;
            FSLib->fsY.ncoor_refine = PETSC_NULL;
        }
    }

    PetscFunctionReturn(0);
}
