// Created by: Jeff Andrews
// Written: Oct 11, 2013
// Last Edited: Oct 6, 2014

/*
This program is designed to take a sample of WDs selected from SDSS
and extract a WD luminosity function. This requires a position and
distance for each WD. The distance is determined photometrically,
by assuming log g = 8.0, and matching to Bergeron's DA/DB models.
*/

#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
#include<sstream>
#include<vector>
#include"LF.h"
#include"matrix.h"

using namespace std;

#define PI acos(-1.0)

int main (int argc, char* argv[]){

  double ra,dec;
  double uu,gg,rr,ii,zz;
  double uu_err,gg_err,rr_err,ii_err,zz_err;
  double ug,gr,ri,iz;
  double pm_ra,pm_dec,pm_tot;
  double pm_ra_err,pm_dec_err;
  double model_temp[60],model_age[60],model_mass[60],model_bol[60];
  double model_u[60],model_g[60],model_r[60],model_i[60],model_z[60];
  double temperature,age,mass,distance,Mbol;
  double l,b,gal_x,gal_y,gal_z;

  string line;
  vector<string> data(20,"0");

  fstream IN;
  ofstream OUT;

  IN.open("./3_rpm.dat");
  OUT.open("./LF.dat");

  create_grids(model_temp,model_age,model_mass,model_bol,model_u,model_g,model_r,model_i,model_z);

  while(getline(IN,line)){

    data.clear();
    split(line, data);

    ra = strtof(data[0]);
    dec = strtof(data[1]);
    uu = strtof(data[2]);
    gg = strtof(data[3]);
    rr = strtof(data[4]);
    ii = strtof(data[5]);
    zz = strtof(data[6]);
    uu_err = strtof(data[7]);
    gg_err = strtof(data[8]);
    rr_err = strtof(data[9]);
    ii_err = strtof(data[10]);
    zz_err = strtof(data[11]);
    pm_ra = strtof(data[12]);
    pm_dec = strtof(data[13]);
    pm_ra_err = strtof(data[14]);
    pm_dec_err = strtof(data[15]);

    calc_wd(model_temp,model_age,model_mass,model_bol,model_u,model_g,model_r,model_i,model_z,uu,gg,rr,ii,zz,&temperature,&mass,&age,&distance,&Mbol);

    equat_to_galac(ra,dec,&l,&b);
    galac_to_cartes(l,b,distance,&gal_x,&gal_y,&gal_z);

    pm_tot = sqrt(pm_ra*pm_ra*cos(DtoR(dec))*cos(DtoR(dec)) + pm_dec*pm_dec);

    OUT << ra << " " << dec << " ";
    OUT << gg << " ";
    OUT << pm_tot << " ";
    OUT << 4.74e-3 * pm_tot * distance * 1.0e3 << " ";   // proper motion in km/s
    OUT << temperature << " ";
    OUT << age << " ";
    OUT << distance << " ";
    OUT << Mbol << " ";
    OUT << RtoD(l) << " ";
    OUT << RtoD(b) << " ";
    OUT << gal_x << " ";
    OUT << gal_y << " ";
    OUT << gal_z << " ";
    OUT << endl;
  }

  


  return 0;
}

void create_grids(double model_temp[60],double model_age[60],double model_mass[60], double model_bol[60], double model_u[60], double model_g[60], double model_r[60], double model_i[60], double model_z[60]){
  int i;

  vector<string> data(30,"0");

  string line,dir;
  string file;
  ifstream IN;

  IN.open("/home/jeff/Research/data/bergeronDA8.dat");

  for(i=0;i<60;i++){
    getline(IN,line);
    split(line,data);

    model_temp[i] = strtof(data[0]);
    model_age[i] = strtof(data[27]);
    model_mass[i] = strtof(data[2]);
    model_bol[i] = strtof(data[3]);
    model_u[i] = strtof(data[13]);
    model_g[i] = strtof(data[14]);
    model_r[i] = strtof(data[15]);
    model_i[i] = strtof(data[16]);
    model_z[i] = strtof(data[17]);
  }

  IN.close();

}

void calc_wd(double model_temp[60],double model_age[60],double model_mass[60],double model_Mbol[60],double model_u[60],double model_g[60],double model_r[60],double model_i[60],double model_z[60],double uu,double gg,double rr,double ii,double zz,double* temperature,double* mass,double* age,double* distance,double* Mbol){
  int i,best_i;
  double best_temp,best_age,best_mass,best_dist,best_Mbol;
  double ug,gr,ri,iz;
  double d_u,d_g,d_r,d_i,d_z;
  double model_ug,model_gr,model_ri,model_iz;
  double temp_u,temp_g,temp_r,temp_i,temp_z;
  double temp_ug,temp_gr,temp_ri,temp_iz;
  double temp, fit, best_fit;

  double x[3];
  double A_age[3],A_m[3],A_bol[3],A_u[3],A_g[3],A_r[3],A_i[3],A_z[3];
  double z_age[3],z_m[3],z_bol[3],z_u[3],z_g[3],z_r[3],z_i[3],z_z[3];

  ug = uu-gg;
  gr = gg-rr;
  ri = rr-ii;
  iz = ii-zz;
  best_fit = 1000.0;
  for(i=0;i<59;i++){
    fit = 0.0;
    model_ug = model_u[i]-model_g[i];
    model_gr = model_g[i]-model_r[i];
    model_ri = model_r[i]-model_i[i];
    model_iz = model_i[i]-model_z[i];

    fit += (model_ug-ug)*(model_ug-ug);
    fit += (model_gr-gr)*(model_gr-gr);
    fit += (model_ri-ri)*(model_ri-ri);
    fit += (model_iz-iz)*(model_iz-iz);

    if (fit < best_fit){
      best_fit = fit;
      best_i = i;
    }
  }

  if (best_i == 0) best_i = 1;
  if (best_i == 59) best_i = 58;
  x[0] = model_temp[best_i-1];
  x[1] = model_temp[best_i];
  x[2] = model_temp[best_i+1];
  z_age[0] = model_age[best_i-1];
  z_age[1] = model_age[best_i];
  z_age[2] = model_age[best_i+1];
  z_m[0] = model_mass[best_i-1];
  z_m[1] = model_mass[best_i];
  z_m[2] = model_mass[best_i+1];
  z_bol[0] = model_Mbol[best_i-1];
  z_bol[1] = model_Mbol[best_i];
  z_bol[2] = model_Mbol[best_i+1];
  z_u[0] = model_u[best_i-1];
  z_u[1] = model_u[best_i];
  z_u[2] = model_u[best_i+1];
  z_g[0] = model_g[best_i-1];
  z_g[1] = model_g[best_i];
  z_g[2] = model_g[best_i+1];
  z_r[0] = model_r[best_i-1];
  z_r[1] = model_r[best_i];
  z_r[2] = model_r[best_i+1];
  z_i[0] = model_i[best_i-1];
  z_i[1] = model_i[best_i];
  z_i[2] = model_i[best_i+1];
  z_z[0] = model_z[best_i-1];
  z_z[1] = model_z[best_i];
  z_z[2] = model_z[best_i+1];

  find_coeff_3(A_age,x,z_age);
  find_coeff_3(A_m,x,z_m);
  find_coeff_3(A_bol,x,z_bol);
  find_coeff_3(A_u,x,z_u);
  find_coeff_3(A_g,x,z_g);
  find_coeff_3(A_r,x,z_r);
  find_coeff_3(A_i,x,z_i);
  find_coeff_3(A_z,x,z_z);

  best_fit = 1000.0;  
  for(i=0;i<1000;i++){
    fit = 0.0;
    temp = (double)i*(x[2]-x[0])/1000.0 + x[0];

    temp_u = interp_3(temp,A_u);
    temp_g = interp_3(temp,A_g);
    temp_r = interp_3(temp,A_r);
    temp_i = interp_3(temp,A_i);
    temp_z = interp_3(temp,A_z);

    temp_ug = temp_u-temp_g;
    temp_gr = temp_g-temp_r;
    temp_ri = temp_r-temp_i;
    temp_iz = temp_i-temp_z;

    fit += (temp_ug-ug)*(temp_ug-ug);
    fit += (temp_gr-gr)*(temp_gr-gr);
    fit += (temp_ri-ri)*(temp_ri-ri);
    fit += (temp_iz-iz)*(temp_iz-iz);
  
    if(fit < best_fit){
      best_fit = fit;
      best_temp = temp;
      best_age = interp_3(temp,A_age);
      best_mass = interp_3(temp,A_m);
      best_Mbol = interp_3(temp,A_bol);

      d_u = pow(10.0,0.2*(uu-temp_u-10.0));
      d_g = pow(10.0,0.2*(gg-temp_g-10.0));
      d_r = pow(10.0,0.2*(rr-temp_r-10.0));
      d_i = pow(10.0,0.2*(ii-temp_i-10.0));
      d_z = pow(10.0,0.2*(zz-temp_z-10.0));

      best_dist = (d_u+d_g+d_r+d_i+d_z)/5.0;
    }

  }

  *temperature = best_temp;
  *age = best_age;
  *mass = best_mass;
  *distance = best_dist;
  *Mbol = best_Mbol;
}

void equat_to_galac(double ra,double dec,double* l,double* b){
  double ra_J,dec_J,ra_B,dec_B;
  double temp1,temp2;
  double Bpos[3],Jpos[3];
  double BtoJ[3][3];
  double JtoB[3][3];

  
  BtoJ[0][0] = 0.999926;
  BtoJ[0][1] = -0.011179;
  BtoJ[0][2] = -0.004859;
  BtoJ[1][0] = 0.011179;
  BtoJ[1][1] = 0.999938;
  BtoJ[1][2] = -0.000027;
  BtoJ[2][0] = 0.004859;
  BtoJ[2][1] = 0.000027;
  BtoJ[2][2] = 0.999988;

  Jpos[0] = cos(DtoR(ra))*cos(DtoR(dec));
  Jpos[1] = sin(DtoR(ra))*cos(DtoR(dec));
  Jpos[2] = sin(DtoR(dec));

  matrix_iden_3(JtoB);
  invert_3(BtoJ,JtoB);

  matrix_mult_3(Bpos,Jpos,JtoB);
  
  ra_B = atan(Bpos[1]/Bpos[0]);
  dec_B = asin(Bpos[2]);

  *b = asin(0.4602*sin(dec_B) - 0.8878*cos(dec_B)*sin(ra_B-4.9262));
  *l = asin((0.4602*cos(dec_B)*sin(ra_B-4.9262) + 0.8878*sin(dec_B))/cos(*b));
  *l = *l + 0.57596;

}

void galac_to_cartes(double l,double b,double distance,double* gal_x,double* gal_y,double* gal_z){
  double dist_sun_to_center = 7.8;
  double dist_sun_off_plane = 0.026;

  *gal_x = -distance*cos(l)*cos(b) + dist_sun_to_center;
  *gal_y = distance*sin(l)*cos(b);
  *gal_z = distance*sin(b) + dist_sun_off_plane;
}

double DtoR(double theta_d){
  return theta_d*PI/180.0;
}

double RtoD(double theta_r){
  return theta_r*180.0/PI;
}

int strtoi(string temp){
  int out;
  istringstream str_in(temp);
  str_in >> out;
  return out;
}

string strtoa(string temp){
  string out;
  istringstream str_in(temp); 
  str_in >> out;
  return out;
}

double strtof(string temp){
  double out;
  istringstream str_in(temp); 
  str_in >> out;
  return out;
}

void split(string line, vector<string>& str_out){
  string tmp;
  size_t spot,spot_s,spot_t;
  str_out.clear();
  while(!line.empty()){
    spot_s = line.find(" ");
    spot_t = line.find("\t");

    if (spot_s == 0 || spot_t == 0){
      line = line.substr(1,line.size());
      continue;
    }

    if(spot_s<1000){
      if(spot_t>1000){
        spot = spot_s;
      }else {
        spot = min(spot_s,spot_t);
      }
    }else if(spot_t<1000){
      spot = spot_t;
    }else {
      spot = line.size();  // For the last element
    }

    tmp = line.substr(0,spot);
    str_out.push_back(tmp);
    tmp.clear();
    line.erase(0,spot+1);
  }
}


