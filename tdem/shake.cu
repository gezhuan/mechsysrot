/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2021 Sergio Galindo                                    *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
*                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/
// Periodic boundary conditions


// MechSys

#include <math.h>
#include <iostream>
#include <fstream>
#include <mechsys/dem/domain.h>
using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
    size_t Nproc = 0.75*omp_get_max_threads();

    String filekey  (argv[1]);
	String filename (filekey+".inp");
	if (!Util::FileExists(filename)) throw new Fatal("File <%s> not found",filename.CStr());
	ifstream infile(filename.CStr());

    double alphaprop = 2.0;
    if (argc>=3) Nproc         = atoi(argv[2]);
    if (argc>=4) alphaprop     = atof(argv[3]);


    String ptype;
    bool   Render = true;
    double L_x       ;
    double L_y       ;
    double L_z       ;
    double nu        ;
    double dx        ;
    double dt        ;
    double Dp        ;
    double R         ;   
    double Tf1       ;
    double Tf2       ;
    double dtOut1    ;
    double dtOut2    ;
    double height    ;
    double rho_s     ;
    double rho_l     ;
    double Mu        ;
    double Muw       ;
    double Beta      ;        //
	double Kn        ;         // Normal stiffness
	double Kt        ;         // Tangential stiffness
	double Gn        ;         // Normal dissipative coefficient
	double ratiots   ;     // Coefficient to modify the time step
	double ratiokn   ;
    double R_random  ;
    double Gnw       ;
    double sizedis    ;
    double MUP       ;
    double Ang       ;
    double fraction  ;
    int    Nump      ;
    double Amp       ;
    double Tper      ;



    
    {
        infile >> ptype;           infile.ignore(200,'\n');
        infile >> Render;          infile.ignore(200,'\n');
        infile >> L_x;             infile.ignore(200,'\n');
        infile >> L_y;             infile.ignore(200,'\n');
        infile >> L_z;             infile.ignore(200,'\n');
        infile >> nu;              infile.ignore(200,'\n');
        infile >> dx;              infile.ignore(200,'\n');
        infile >> dt;              infile.ignore(200,'\n');
        infile >> Dp;              infile.ignore(200,'\n');
        infile >> R;               infile.ignore(200,'\n');
        infile >> Tf1;             infile.ignore(200,'\n');
        infile >> Tf2;             infile.ignore(200,'\n');
        infile >> dtOut1;          infile.ignore(200,'\n');
        infile >> dtOut2;          infile.ignore(200,'\n');
        infile >> height;          infile.ignore(200,'\n');
        infile >> rho_s;           infile.ignore(200,'\n');
        infile >> rho_l;           infile.ignore(200,'\n');
        infile >> Mu;              infile.ignore(200,'\n');
        infile >> Muw;             infile.ignore(200,'\n');
        infile >> Beta;            infile.ignore(200,'\n');
        infile >> Kn;              infile.ignore(200,'\n');
        infile >> Gn;              infile.ignore(200,'\n');
        infile >> ratiots;         infile.ignore(200,'\n');
        infile >> ratiokn;         infile.ignore(200,'\n');
        infile >> R_random;        infile.ignore(200,'\n');
        infile >> Gnw;             infile.ignore(200,'\n');
        infile >> sizedis;          infile.ignore(200,'\n');
        infile >> MUP;             infile.ignore(200,'\n');
        infile >> Ang;             infile.ignore(200,'\n');
        infile >> fraction;             infile.ignore(200,'\n');
        infile >> Nump;             infile.ignore(200,'\n');
        infile >> Amp;             infile.ignore(200,'\n');
        infile >> Tper;             infile.ignore(200,'\n');

    }

    
    Kt=ratiokn*Kn;
  

    double acc=Dp * (rho_s - rho_l) / rho_s;
    double ang=Ang*M_PI/180.0;
    cout <<"acc"<< acc<<endl;
    cout <<"L_z"<< L_z<<endl;
    double Sphere_size;
    double dtdem;
    Vec3_t Xmin0;
    Vec3_t Xmax0;
    Vec3_t Center;
    double R_base=9.2;
    Array<int> delpar0;
    Vec3_t axis0(OrthoSys::e0); // rotation of face
    Vec3_t axis1(OrthoSys::e1); // rotation of face
    double thichness=4*R;
    


    DEM::Domain d2;
    
    Center = Vec3_t(0.0,0.0,0.0);
    
    
    d2.AddSphere(-1,Vec3_t(Center(0)-2.0*R,Center(1),Center(2)+2.05*R),R,rho_s);
    d2.AddSphere(-3,Vec3_t(Center(0)-2.0*R,Center(1),Center(2)),R,rho_s);

    d2.AddSphere(-2,Vec3_t(Center(0)+2.0*R,Center(1),Center(2)+2.05*R),R,rho_s);

    d2.AddSphere(-4,Vec3_t(Center(0)+2.0*R,Center(1),Center(2)),R,rho_s);

    double angle = M_PI*Tper;
    Vec3_t axis(1.0,0.0,0.0);
    Quaternion_t q;
    NormalizeRotation (angle,axis,q);
    d2.GetParticle(-1)->Q = q;
    d2.GetParticle(-4)->Q = q;

    d2.GetParticle(-1)->Props.fac = 1.0e-4;
    d2.GetParticle(-2)->Props.fac = 1.0e-4;
    d2.GetParticle(-3)->Props.fac = 1.0e-4;
    d2.GetParticle(-4)->Props.fac = 1.0;
    

    for (size_t np=0;np<d2.Particles.Size();np++)
    {
        if (d2.Particles[np]->Tag >= -2)
        {
            d2.Particles[np]->v = Vec3_t(0.0,0.0,-2.0);
            d2.Particles[np]->Props.Kn = Kn; // normal stiffness
            d2.Particles[np]->Props.Kt = Kt; // trangential stiffness
            d2.Particles[np]->Props.Gn = Gn; // restitution coefficient
            d2.Particles[np]->Props.Mu = MUP; // frictional coefficient
        }
        else
        {
            d2.Particles[np]->Props.Kn = Kn; // normal stiffness
            d2.Particles[np]->Props.Kt = Kt; // trangential stiffness
            d2.Particles[np]->Props.Gn = Gn; // restitution coefficient
            d2.Particles[np]->Props.Mu = MUP; // frictional coefficient
        }

    }

    dtdem = ratiots*d2.CriticalDt(); //Calculating time step
    d2.Alpha = R/4.0; //Verlet distance
    d2.Solve(/*tf*/Tf2, dtdem, /*dtOut*/dtOut2, NULL, NULL, "shake", 2, Nproc);
    d2.Save("Stage_P2");
    

}
MECHSYS_CATCH

