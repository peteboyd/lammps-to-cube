#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define KNRM  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"
#define ANGS2BOHR 1.88973

typedef enum { false, true } bool;

typedef struct cbf{
    int natoms, ntypes, npoints, ngx, ngy, ngz;
    double ox, oy, oz;
    double dgx, dgy, dgz;
    double *masses;
    double *xxx;
    double *yyy;
    double *zzz;
    double *charge;
    double *cella;
    double *cellb;
    double *cellc;
    double *icella;
    double *icellb;
    double *icellc;
    double *grid;
    int *types;
    int *alltypes;

} cube_file;

void deconstruct_cube(cube_file *cb){
    free(cb->xxx);
    free(cb->yyy);
    free(cb->zzz);
    free(cb->grid);
    free(cb->masses);
    free(cb->charge);
    free(cb->cella);
    free(cb->cellb);
    free(cb->cellc);
    free(cb->icella);
    free(cb->icellb);
    free(cb->icellc);
    free(cb->types);
    free(cb->alltypes);   
}
void parse_datafile(FILE* data_file, cube_file *cb);
void parse_histfile(FILE* hist_file, cube_file *cb);
void write_cubefile(FILE* cub_file, cube_file *cb);
void deconstruct_cube(cube_file *cb);
int parse_line(char* line, char*** res);
int to_atomic_number(double);

void construct_cell_from_lammps(double, double, double, 
                                double, double, double, 
                                double, double, double, cube_file*);

double* reduced_coord(int, cube_file*);
double* cross(double*, double*);
double dot(double*, double*);

int main(int argc, char**argv)
{
    FILE *df;
    FILE *hf;
    FILE *cf;
    int i;
    if (argc != 3){
        printf("usage: %s [data. file] [lammps histogram file]\n", argv[0]);
        exit(0);
    }
    char const* const dataFileName = argv[1];
    char const* const histFileName = argv[2];
    cube_file cb;

    cb.ngx = 0;
    cb.ngy = 0;
    cb.ngz = 0;
    cb.dgx = 0.0;
    cb.dgy = 0.0;
    cb.dgz = 0.0;

    df = fopen(dataFileName, "r");
    hf = fopen(histFileName, "r");
    cf = fopen("out.cube", "w");
    parse_datafile(df, &cb);
    parse_histfile(hf, &cb);

    write_cubefile(cf, &cb);

    printf("%sCube file written to %s\n%s",KRED,"out.cube",KNRM);
    fclose(df);
    fclose(hf);
    fclose(cf);
    
    deconstruct_cube(&cb);
    return 0;
}

int to_atomic_number(double mass){
    /*
    1  12.0107 # C_R
    2   1.0079 # H_
    3  14.0067 # N_R
    4  14.0067 # N_3
    5  15.9994 # O_R
    6  65.3800 # Zn3+2
    7  12.0107 # C_3
    */
    double tol = 0.2;
    //printf("%f - 14.0 = %f \n", mass, fabs(mass-14.0));
    if (fabs(mass-14.0) < tol){
        return 7;
    }
    else if (fabs(mass - 12.0) < tol){
        return 6;
    }
    else if (fabs(mass - 65.38) < tol){
        return 30;
    }
    else if (fabs(mass - 1.00) < tol){
        return 1;
    }
    else if (fabs(mass - 15.9994) < tol){
        return 8;
    }
    return 0;
}

int parse_line(char* line, char*** res){
    int n_spaces=0;
    int i;
    char* p;
    p = strtok(line, " ");
    while (p){
        *res = realloc(*res, sizeof(char*) * ++n_spaces);
        if (res == NULL)
            exit(-1);
        //strcpy(res[n_spaces - 1], p);
        (*res)[n_spaces-1] = p;
        //Ignore comments.
        //printf("%s n_spaces = %i, p = %s %s\n",KRED,n_spaces, p,KNRM);
        if (strstr(p,"#") != NULL){
            n_spaces--;
            p = 0;
        }
        else{
            p = strtok(NULL, " ");
        }
    }
    *res = realloc(*res, sizeof(char*) *(n_spaces+1));
    (*res)[n_spaces] = 0;
    //printf("%s %s %s\n",KMAG, res[0], KNRM);
    return n_spaces;
}

double *reduced_coord(int i, cube_file *cb){
    double *ans;
    ans = malloc(sizeof(double)*3);
    ans[0] = cb->xxx[i]*cb->icella[0] + cb->yyy[i]*cb->icellb[0] + cb->zzz[i]*cb->icellc[0];
    ans[1] = cb->xxx[i]*cb->icella[1] + cb->yyy[i]*cb->icellb[1] + cb->zzz[i]*cb->icellc[1];
    ans[2] = cb->xxx[i]*cb->icella[2] + cb->yyy[i]*cb->icellb[2] + cb->zzz[i]*cb->icellc[2];
    return ans;
}

double dot(double *avec, double *bvec){
    return (avec[0]*bvec[0] + avec[1]*bvec[1] + avec[2]*bvec[2]);
}

double *cross(double *avec, double *bvec){
    double* ans = malloc(sizeof(double)*3);

    ans[0] = avec[1]*bvec[2] - avec[2]*bvec[1];
    ans[1] = avec[2]*bvec[0] - avec[0]*bvec[2];
    ans[2] = avec[0]*bvec[1] - avec[1]*bvec[0];
    return ans;

}
void construct_cell_from_lammps(double xlo, double xhi,
                                double ylo, double yhi,
                                double zlo, double zhi,
                                double xy, double xz, double yz,
                                cube_file* cb){

    double vol;
    cb->cella = malloc(sizeof(double)*3);
    cb->cellb = malloc(sizeof(double)*3);
    cb->cellc = malloc(sizeof(double)*3);

    const double avec[3] = {xhi - xlo, 0.0, 0.0};
    const double bvec[3] = {xy, yhi-ylo, 0.0};
    const double cvec[3] = {xz, yz, zhi-zlo};

    memcpy(cb->cella, avec, sizeof(avec));
    memcpy(cb->cellb, bvec, sizeof(bvec));
    memcpy(cb->cellc, cvec, sizeof(cvec));


    //inverse cell
    vol = abs(dot(cross(cb->cellc, cb->cella), cb->cellb));
    cb->icella = malloc(sizeof(double)*3);
    cb->icellb = malloc(sizeof(double)*3);
    cb->icellc = malloc(sizeof(double)*3);

    cb->icella[0] = (bvec[1]*cvec[2] - bvec[2]*cvec[1])/vol;
    cb->icella[1] = (avec[2]*cvec[1] - avec[1]*cvec[2])/vol;
    cb->icella[2] = (avec[1]*bvec[2] - avec[2]*bvec[1])/vol;

    cb->icellb[0] = (bvec[2]*cvec[0] - bvec[0]*cvec[2])/vol;
    cb->icellb[1] = (avec[0]*cvec[2] - avec[2]*cvec[0])/vol;
    cb->icellb[2] = (avec[2]*bvec[0] - avec[0]*bvec[2])/vol;

    cb->icellc[0] = (bvec[0]*cvec[1] - bvec[1]*cvec[0])/vol;
    cb->icellc[1] = (avec[1]*cvec[0] - avec[0]*cvec[1])/vol;
    cb->icellc[2] = (avec[0]*bvec[1] - avec[1]*bvec[0])/vol;


    printf("%s Volume of system = %.2f Angstroms^3 %s\n",KRED,vol,KNRM);
}
void parse_histfile(FILE* hist_file, cube_file* cb){
    bool gridread=false;
    bool gxdone=false;
    bool gydone=false;
    bool gzdone=false;
    char line[500], cpline[500];
    char **spline = NULL;
    int i=0, n_spaces;
    int gzcount=0;
    int gycount=0;
    int gxcount=0;
    double pgx, pgy, pgz, gx,gy,gz;

    while(fgets(line, sizeof(line), hist_file)){
        strcpy(cpline, line);
        n_spaces = parse_line(cpline, &spline);
        if((gridread)&&(i<cb->npoints)){
            // n_spaces will indicate how many values are on each line of data.
            // This will have to be adjusted by hand depending on the output of LAMMPS.
            if(n_spaces == 6){
                i++;
                cb->grid[i-1] = atof(spline[4]);
                gx = atof(spline[1]);
                gy = atof(spline[2]);
                gz = atof(spline[3]);
                if(i==1){
                    cb->ox = gx; 
                    cb->oy = gy;
                    cb->oz = gz;
                }
                else{
                    if((gx-pgx != 0.0)&&(cb->dgx == 0.0)){
                        gxcount++;
                        cb->dgx = gx - pgx;
                    }
                    else if(gx -pgx != 0.0){
                        gxcount++;
                    }
                    if((gx - pgx != 0.0)&&(gx - cb->ox == 0.0)&&(cb->dgx != 0.0)){
                        cb->ngx = gxcount;
                    }
                    if(!gydone){
                        if((gy-pgy != 0.0)&&(cb->dgy == 0.0)){
                            gycount++;
                            cb->dgy = gy - pgy;
                        }
                        else if(gy-pgy !=0.0){
                            gycount++;
                        }
                        if((gy-pgy != 0.0)&&(gy-cb->oy == 0.0)&&(cb->dgy != 0.0)){
                            //printf("%.2f - %.2f = %.2f\n", gy, pgy, gy-pgy);
                            //printf("step %i, dgy =  %.2f \n", i, cb->dgy);
                            cb->ngy = gycount;
                            gydone=true;
                        }
                    }
                    if(!gzdone){
                        if((gz-pgz != 0.0)&&(cb->dgz == 0.0)){
                            gzcount++;
                            cb->dgz = gz - pgz;
                        }
                        else if(gz-pgz !=0.0){
                            gzcount++;
                        }
                        if((gz-pgz != 0.0)&&(gz - cb->oz == 0.0)&&(cb->dgz != 0.0)){
                            cb->ngz = gzcount;
                            gzdone=true;
                        }
                    }
                }
                pgx = gx;
                pgy = gy;
                pgz = gz;
            }
        }
        if((gridread)&&(i >= cb->npoints)){
            cb->ngx = gxcount + 1;
            gycount = 0;
            gzcount = 0;
            gridread=false;
        }
        if(n_spaces == 3){
            cb->npoints = atoi(spline[1]);
            cb->grid = malloc(sizeof(double)*cb->npoints);
            gridread=true;
        }
    }

}

void parse_datafile(FILE* data_file, cube_file* cb){
    /*
     * assumes a monoclinic cell
           0.000000  30.409200 xlo xhi
           0.000000  30.633400 ylo yhi
           0.000000  43.114053 zlo zhi
           0.000000   0.000000   0.112496 xy xz yz
    */
    bool atread=false, typeread=false;
    double xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz;
    char line[500],cpline[500];
    char **spline=NULL;
    int i=0, n_spaces;

    while (fgets(line, sizeof(line), data_file)){
        strcpy(cpline, line);
        n_spaces = parse_line(cpline, &spline);
        
        //read in atoms    
        if((atread)&&(i<cb->natoms)){
            if(n_spaces == 7){
                i++;
                cb->alltypes[i-1] = atoi(spline[2]);
                cb->charge[i-1] = atof(spline[3]);
                cb->xxx[i-1] = atof(spline[4]);
                cb->yyy[i-1] = atof(spline[5]);
                cb->zzz[i-1] = atof(spline[6]);

            }
        }
        else if((atread) && (i>=cb->natoms)){
            i=0;
            atread=false;
        }
        //read in atom types
        if((typeread) && (i < cb->ntypes)){
            if(n_spaces == 2){
                i++;
                cb->types[i-1] = atoi(spline[0]);
                cb->masses[i-1] = atof(spline[1]);
            }
        }
        else if((typeread) && (i>=cb->ntypes)){
            i=0;
            typeread=false;
        }
        if(strstr(line, "atoms") != NULL){
            cb->natoms = atoi(spline[0]);
            cb->xxx = malloc(sizeof(double)*cb->natoms+1);
            cb->yyy = malloc(sizeof(double)*cb->natoms+1);
            cb->zzz = malloc(sizeof(double)*cb->natoms+1);
            cb->charge = malloc(sizeof(double)*cb->natoms+1);
            cb->alltypes = malloc(sizeof(int)*cb->natoms+1);

        }
        else if(strstr(line, "atom types") != NULL){
            cb->ntypes=atoi(spline[0]);
            cb->types = malloc(sizeof(int) * cb->ntypes+1);
            cb->masses = malloc(sizeof(double) * cb->ntypes+1);
        }
        else if(strstr(line, "xlo") != NULL){
            xlo=atof(spline[0]);
            xhi=atof(spline[1]);
        }
        else if(strstr(line, "ylo") != NULL){
            ylo=atof(spline[0]);
            yhi=atof(spline[1]);

        }
        else if(strstr(line, "zlo") != NULL){
            zlo=atof(spline[0]);
            zhi=atof(spline[1]);

        }
        else if(strstr(line, " xy ") != NULL){
            xy = atof(spline[0]);
            xz = atof(spline[1]);
            yz = atof(spline[2]);
        }
        else if(strstr(line, "Atoms") != NULL){
            atread=true;
        }
        else if(strstr(line, "Masses") != NULL){
            typeread=true;
        }
    }
    free(spline);
    construct_cell_from_lammps(xlo, xhi, ylo, 
                               yhi, zlo, zhi,
                               xy, xz, yz, cb);

}

void write_cubefile(FILE* cb_file, cube_file* cb){
    double *redcrd;
    double mass; 
    int i, j, atmnbr, type;
    printf("%s Number of atom types: %i %s\n",KCYN, cb->ntypes, KNRM);
    fprintf(cb_file,"---------------------LAMMPS CUBE FILE---------------------\n");
    fprintf(cb_file,"-------cube file created from the Lammps data files-------\n");
    printf("%s Number of atoms: %i %s\n",KCYN, cb->natoms, KNRM);
    fprintf(cb_file,"%6i %12.6f %12.6f %12.6f\n",cb->natoms, cb->ox, cb->oy, cb->oz);
    printf("%s Number of grid points: %i %s\n",KCYN, cb->npoints, KNRM);
    printf("%s Number of grid points in x direction: %i %s\n",KGRN, cb->ngx, KNRM);
    fprintf(cb_file, "%6i %12.6f %12.6f %12.6f \n",cb->ngx,
                                                   cb->cella[0]/(double)cb->ngx*ANGS2BOHR, 
                                                   cb->cella[1]/(double)cb->ngx*ANGS2BOHR, 
                                                   cb->cella[2]/(double)cb->ngx*ANGS2BOHR);
    printf("%s Number of grid points in y direction: %i %s\n",KYEL, cb->ngy, KNRM);
    fprintf(cb_file,"%6i %12.6f %12.6f %12.6f \n",cb->ngy,
                                                  cb->cellb[0]/(double)cb->ngy*ANGS2BOHR, 
                                                  cb->cellb[1]/(double)cb->ngy*ANGS2BOHR, 
                                                  cb->cellb[2]/(double)cb->ngy*ANGS2BOHR);
    printf("%s Number of grid points in z direction: %i %s\n",KRED, cb->ngz, KNRM);
    fprintf(cb_file, "%6i %12.6f %12.6f %12.6f \n",cb->ngz, 
                                                   cb->cellc[0]/(double)cb->ngz*ANGS2BOHR, 
                                                   cb->cellc[1]/(double)cb->ngz*ANGS2BOHR, 
                                                   cb->cellc[2]/(double)cb->ngz*ANGS2BOHR);
    for(i= 0; i<cb->ntypes; i++){
        printf("%s Mass of type %i = %.4f \n%s", KCYN, cb->types[i], cb->masses[i], KNRM);
    }
    for(i = 0; i<cb->natoms; i++){
        type = cb->alltypes[i];
        mass = cb->masses[type-1];
        atmnbr = to_atomic_number(mass); 
        //redcrd = reduced_coord(i, cb);
        fprintf(cb_file,"%6i %12.6f %12.6f %12.6f %12.6f\n", atmnbr, mass, cb->xxx[i]*ANGS2BOHR, 
                                                                           cb->yyy[i]*ANGS2BOHR,
                                                                           cb->zzz[i]*ANGS2BOHR);
        //free(redcrd);
    }
    //Format the grid in columns of 6
    //c goes first (z)
    i=0;
    while(i < cb->npoints){
       for(j=0; j<6; j++){
           if(j==0){
               fprintf(cb_file,"%.6E",cb->grid[i]);
           }
           else{
               fprintf(cb_file,"  %.6E",cb->grid[i]);
           }
           i++;
           if((i % cb->ngz) == 0 ){
               break;
           }
       }
       fprintf(cb_file,"\n");
    }
}
