# M3TET_SPH - 3D Finite Element Convection Code in Spherical Geometry

### Authors:
Jorge M. Taramón (jorge.taramongomez.2014@alumni.rhul.ac.uk), Jason P. Morgan, Jörg Hasenclever

### Description:
- M3TET_SPH solves for thermo-mechanical viscous flow evolution in spherical geometries
- Setup is defined in `model_parameters_sph_south_atlantic_plume.m` so it can simulate the influence of a plume during the South Atlantic rifting
- Velocity BCs are free-slip along the core-mantle boundary and prescribed plate motion on the top surface using plate kinematic reconstructions (Gurnis et al., 2012)
- The unstructured finite element spherical mesh is generated using [springmesh_3d_spherical_shell](https://github.com/JorgeTaramon/Mesh_Generator/tree/master/springmesh_3d_spherical_shell) ([Taramón et al., 2019](https://doi.org/10.1016/j.cageo.2019.104324))
- The Double-Jacobian approach [Morgan et al., 2019](https://onlinelibrary.wiley.com/doi/abs/10.1002/fld.4799) is used to improve the effieciency when solving solving problems in a spherical geometry as well as to speed-up the particle search routines in curved-edge elements

### Running this code for the first time:
1. Download "SuiteSparse" from:
    http://faculty.cse.tamu.edu/davis/SuiteSparse/SuiteSparse-4.0.2.tar.gz
2. Install it by running SuiteSparse_install.m
3. Download mutils-0.4-2 from:
    https://sourceforge.net/projects/milamin/files/mutils-0.4-2.zip
4. Install it by running install.m (follow the instructions that will appear in the Command Window)
5. In order to set the path for mutils
    - 5.1. Write pctconfig in the Command Window
    - 5.2. Open `SETUP_TEST/addpaths_mutils.m`
    - 5.3. Create a new case for your computer:
        - Enter the hostname you obtained in 5.1 as a new ‘case’, e.g., case 'fpdc462'
        - Enter the path for the folder mutils in ‘path2mutils’
6. In order to set the OUTPUT folder in your computer:
    - 6.1. Open `SETUP_TEST/data_storage_3d.m`
    - 6.2. Create a new case for your computer:
        - Enter the hostname you obtained in 5.1 as a new ‘case’, e.g., case 'fpdc462'
        - Enter the path to save the output data in ‘path2data’
7. Run `M3TET_SPH_SOUTH_ATLANTIC_PLUME.m`
