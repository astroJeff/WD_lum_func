#include<string>

using namespace std;

int strtoi(string temp);
string strtoa(string temp);
double strtof(string temp);
void split(string line, vector<string>& str_out); 
void calc_wd(double model_temp[60],double model_age[60],double model_mass[60],double model_Mbol[60],double model_u[60],double model_g[60],double model_r[60],double model_i[60],double model_z[60],double uu,double gg,double rr,double ii,double zz,double* temperature,double* mass,double* age,double* distance,double* Mbol);
void create_grids(double model_temp[60],double model_age[60],double model_mass[60],double model_bol[60],double model_u[60],double model_g[60],double model_r[60],double model_i[60],double model_z[60]);
void equat_to_galac(double ra,double dec,double* l,double* b);
void galac_to_cartes(double l,double b, double distance,double* gal_x,double* gal_y,double* gal_z);
double DtoR(double theta_d);
double RtoD(double theta_r);


