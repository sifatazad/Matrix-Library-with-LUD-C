#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define Abs(a) ((a)>=0?(a):-(a))

static int total_row(FILE *fp) {
    int n=0;
    char line[5000];
    while (fgets(line,5000,fp)!=NULL) {
        n++;
    }
    rewind(fp);
    return n;
}

static void lu_decomposition(FILE *fp, FILE *fo, double *soln) {
    int n= total_row(fp);

    double **matrix= malloc(sizeof(double*)*n);
    if (!matrix) {
        return;
    }
    for (int i=0; i<n; i++) {
        matrix[i]=malloc(sizeof(double)*(n+1));
        if (!matrix[i]) {
            free(matrix);
            return;
        }
    }

    for (int i=0; i<n; i++) {
        for (int j=0; j<n+1; j++) {
            fscanf(fp,"%lf\t",&matrix[i][j]);
        }
    }
    rewind(fp);

    for (int step=0; step<n-1; step++) {
        double col_max=Abs(matrix[step][step]);

        for (int i=step+1; i<n; i++) {
            if (Abs(matrix[i][step])>Abs(col_max)) {
                col_max=Abs(matrix[i][step]);
            }
        }

        if (col_max<1e-10) {
            printf("Singular Matrix. No unique Solution!\n");
            for (int ii=0; ii<n; ii++) {
                free(matrix[ii]);
            }
            free(matrix);
            return;
        }

        if (Abs(matrix[step][step])!=col_max) {
            for (int i=step+1; i<n; i++) {
                if (col_max==Abs(matrix[i][step])) {
                    double *temp= malloc(sizeof(double)*(n+1));
                    if (!temp) {
                        for (int ii=0; ii<n; ii++) {
                            free(matrix[ii]);
                        }
                        free(matrix);
                        return;
                    }

                    for (int it=0; it<n+1; it++) {
                        temp[it]=matrix[i][it];
                    }
                    for (int it=0; it<n+1; it++) {
                        matrix[i][it]=matrix[step][it];
                    }
                    for (int it=0; it<n+1; it++) {
                        matrix[step][it]=temp[it];
                    }
                    // pivoting done. now gaussian kora need
                    free(temp);
                    break;
                }
            }
        }

        // gaussian elimination here
        for (int i=step+1; i<n; i++) {
            double factor= matrix[i][step]/matrix[step][step];
            for(int j=step; j<n; j++) {
                if (j==step) {
                    matrix[i][j]= factor;
                }
                else {
                    matrix[i][j]+= -factor*matrix[step][j];
                }
            }
        }
    }

    if (Abs(matrix[n-1][n-1])<1e-10) {
        printf("Singular Matrix!\n");
        return;
    }

    for (int i=0; i<n; i++) {
        for (int j=0; j<n+1; j++) {
            fprintf(fo,"%.6lf\t",matrix[i][j]);
        }
        if (i!=n-1)
            fprintf(fo,"\n");
    }

    double *y= malloc(sizeof(double)*n);
    if (!y) {
        for (int i=0; i<n; i++) {
            free(matrix[i]);
        }
        free(matrix);
        return;
    }
    for (int i=0; i<n; i++) {
        y[i]=0;
    }

    for (int i=0; i<n; i++) {
        double s=0;
        for (int j=0; j<i; j++) {
            s+=matrix[i][j]*y[j];
        }
        y[i]=matrix[i][n]-s;
    }

    for (int i=0; i<n; i++) {
        soln[i]=0;
    }

    for (int i=n-1; i>=0; i--) {
        double s=0;
        for (int j=n-1; j>=i+1; j--) {
            s+=matrix[i][j]*soln[j];
        }
        soln[i]=(y[i]-s)/matrix[i][i];
    }

    // free matrix
    for (int i=0; i<n; i++) {
        free(matrix[i]);
    }
    free(matrix);
    free(y);
}

static void lu_decomposition_without_file(double **matrix, double *soln, const int n) {
    for (int step=0; step<n-1; step++) {
        double col_max=Abs(matrix[step][step]);

        for (int i=step+1; i<n; i++) {
            if (Abs(matrix[i][step])>Abs(col_max)) {
                col_max=Abs(matrix[i][step]);
            }
        }

        if (col_max<1e-10) {
            printf("Singular Matrix. No unique Solution!\n");
            return;
        }

        if (Abs(matrix[step][step])!=col_max) {
            for (int i=step+1; i<n; i++) {
                if (col_max==Abs(matrix[i][step])) {
                    double *temp= malloc(sizeof(double)*(n+1));
                    if (!temp) {
                        return;
                    }

                    for (int it=0; it<n+1; it++) {
                        temp[it]=matrix[i][it];
                    }
                    for (int it=0; it<n+1; it++) {
                        matrix[i][it]=matrix[step][it];
                    }
                    for (int it=0; it<n+1; it++) {
                        matrix[step][it]=temp[it];
                    }
                    // pivoting done. now gaussian kora need
                    free(temp);
                    break;
                }
            }
        }

        // gaussian elimination here
        for (int i=step+1; i<n; i++) {
            double factor= matrix[i][step]/matrix[step][step];
            for(int j=step; j<n; j++) {
                if (j==step) {
                    matrix[i][j]= factor;
                }
                else {
                    matrix[i][j]+= -factor*matrix[step][j];
                }
            }
        }
    }

    if (Abs(matrix[n-1][n-1])<1e-10) {
        printf("Singular Matrix!\n");
        return;
    }

    double *y= malloc(sizeof(double)*n);
    if (!y) {
        return;
    }
    for (int i=0; i<n; i++) {
        y[i]=0;
    }

    for (int i=0; i<n; i++) {
        double s=0;
        for (int j=0; j<i; j++) {
            s+=matrix[i][j]*y[j];
        }
        y[i]=matrix[i][n]-s;
    }

    for (int i=0; i<n; i++) {
        soln[i]=0;
    }

    for (int i=n-1; i>=0; i--) {
        double s=0;
        for (int j=n-1; j>=i+1; j--) {
            s+=matrix[i][j]*soln[j];
        }
        soln[i]=(y[i]-s)/matrix[i][i];
    }

    free(y);
}

static void inverse_matrix(FILE *fp, FILE *fo) {
    int n= total_row(fp);
    double **matrix= malloc(sizeof(double*)*n);
    if (!matrix) {
        return;
    }
    for (int i=0; i<n; i++) {
        matrix[i]=malloc(sizeof(double)*n);
        if (!matrix[i]) {
            free(matrix);
            return;
        }
    }

    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            fscanf(fp,"%lf\t",&matrix[i][j]);
        }
    }
    rewind(fp);

    double **soln= malloc(sizeof(double*)*n);
    if (!soln) {
        for (int ii=0; ii<n; ii++) {
            free(matrix[ii]);
        }
        free(matrix);
        return;
    }
    for (int i=0; i<n; i++) {
        soln[i]=malloc(sizeof(double)*n);
        if (!soln[i]) {
            for (int ii=0; ii<n; ii++) {
                free(matrix[ii]);
            }
            free(matrix);
            free(soln);
            return;
        }
    }

    double **aug= malloc(sizeof(double*)*n);
    if (!aug) {
        for (int ii=0; ii<n; ii++) {
            free(matrix[ii]);
        }
        free(matrix);
        for (int i=0; i<n; i++) {
            free(soln[i]);
        }
        free(soln);
        return;
    }
    for (int i=0; i<n; i++) {
        aug[i]=malloc(sizeof(double)*(n+1));
        if (!aug[i]) {
            for (int ii=0; ii<n; ii++) {
                free(matrix[ii]);
            }
            free(matrix);
            for (int ii=0; ii<n; ii++) {
                free(soln[ii]);
            }
            free(soln);
            free(aug);
            return;
        }
    }

    for (int step=0; step<n; step++) {
        double *temp_soln= malloc(sizeof(double)*n);
        if (!temp_soln) {
            for (int ii=0; ii<n; ii++) {
                free(matrix[ii]);
            }
            free(matrix);
            for (int i=0; i<n; i++) {
                free(soln[i]);
            }
            free(soln);
            for (int i=0; i<n; i++) {
                free(aug[i]);
            }
            free(aug);
            return;
        }

        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                aug[i][j]=matrix[i][j];
            }
        }
        for (int i=0; i<n; i++) {
            if (i==step) {
                aug[i][n]=1;
            }
            else {
                aug[i][n]=0;
            }
        }

        lu_decomposition_without_file(aug,temp_soln,n);
        for (int i=0; i<n; i++) {
            soln[i][step]=temp_soln[i];
        }
        free(temp_soln);
    }

    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            fprintf(fo,"%.6lf\t",soln[i][j]);
        }
        if (i!=n-1)
            fprintf(fo,"\n");
    }

    for (int i=0; i<n; i++) {
        free(matrix[i]);
    }
    free(matrix);
    for (int i=0; i<n; i++) {
        free(soln[i]);
    }
    free(soln);
    for (int i=0; i<n; i++) {
        free(aug[i]);
    }
    free(aug);
}

static void inverse_matrix_without_file(double **matrix, double **soln, const int n) {
    double **aug= malloc(sizeof(double*)*n);
    if (!aug) {
        return;
    }
    for (int i=0; i<n; i++) {
        aug[i]=malloc(sizeof(double)*(n+1));
        if (!aug[i]) {
            free(aug);
            return;
        }
    }

    for (int step=0; step<n; step++) {
        double *temp_soln= malloc(sizeof(double)*n);
        if (!temp_soln) {
            for (int i=0; i<n; i++) {
                free(aug[i]);
            }
            free(aug);
            return;
        }

        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                aug[i][j]=matrix[i][j];
            }
        }
        for (int i=0; i<n; i++) {
            if (i==step) {
                aug[i][n]=1;
            }
            else {
                aug[i][n]=0;
            }
        }

        lu_decomposition_without_file(aug,temp_soln,n);
        for (int i=0; i<n; i++) {
            soln[i][step]=temp_soln[i];
        }
        free(temp_soln);
    }

    for (int i=0; i<n; i++) {
        free(aug[i]);
    }
    free(aug);
}

static double determinant(FILE *fp) {
    int n= total_row(fp);

    double **matrix= malloc(sizeof(double*)*n);
    if (!matrix) {
        return NAN;
    }
    for (int i=0; i<n; i++) {
        matrix[i]=malloc(sizeof(double)*n);
        if (!matrix[i]) {
            free(matrix);
            return NAN;
        }
    }

    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            fscanf(fp,"%lf\t",&matrix[i][j]);
        }
    }
    rewind(fp);

    int swap=0;
    for (int step=0; step<n-1; step++) {
        double col_max=Abs(matrix[step][step]);

        for (int i=step+1; i<n; i++) {
            if (Abs(matrix[i][step])>Abs(col_max)) {
                col_max=Abs(matrix[i][step]);
            }
        }

        if (col_max<1e-10) {
            for (int ii=0; ii<n; ii++) {
                free(matrix[ii]);
            }
            free(matrix);
            return 0.0;
        }

        if (Abs(matrix[step][step])!=col_max) {
            swap++;
            for (int i=step+1; i<n; i++) {
                if (col_max==Abs(matrix[i][step])) {
                    double *temp= malloc(sizeof(double)*(n+1));
                    if (!temp) {
                        for (int ii=0; ii<n; ii++) {
                            free(matrix[ii]);
                        }
                        free(matrix);
                        return NAN;
                    }

                    for (int it=0; it<n; it++) {
                        temp[it]=matrix[i][it];
                    }
                    for (int it=0; it<n; it++) {
                        matrix[i][it]=matrix[step][it];
                    }
                    for (int it=0; it<n; it++) {
                        matrix[step][it]=temp[it];
                    }
                    // pivoting done. now gaussian kora need
                    free(temp);
                    break;
                }
            }
        }

        // gaussian elimination here
        for (int i=step+1; i<n; i++) {
            double factor= matrix[i][step]/matrix[step][step];
            for(int j=step; j<n; j++) {
                if (j==step) {
                    matrix[i][j]= factor;
                }
                else {
                    matrix[i][j]+= -factor*matrix[step][j];
                }
            }
        }
    }

    if (Abs(matrix[n-1][n-1])<1e-10) {
        for (int ii=0; ii<n; ii++) {
            free(matrix[ii]);
        }
        free(matrix);
        return 0.0;
    }

    double det=1;
    for (int i=0; i<n; i++) {
        det*=matrix[i][i];
    }

    // free matrix
    for (int i=0; i<n; i++) {
        free(matrix[i]);
    }
    free(matrix);
    return det*pow(-1,swap);
}

static double determinant_without_file(double **matrix, const int n) {
    int swap=0;
    for (int step=0; step<n-1; step++) {
        double col_max=Abs(matrix[step][step]);

        for (int i=step+1; i<n; i++) {
            if (Abs(matrix[i][step])>Abs(col_max)) {
                col_max=Abs(matrix[i][step]);
            }
        }

        if (col_max<1e-10) {
            return 0.0;
        }

        if (Abs(matrix[step][step])!=col_max) {
            swap++;
            for (int i=step+1; i<n; i++) {
                if (col_max==Abs(matrix[i][step])) {
                    double *temp= malloc(sizeof(double)*(n+1));
                    if (!temp) {
                        return NAN;
                    }

                    for (int it=0; it<n; it++) {
                        temp[it]=matrix[i][it];
                    }
                    for (int it=0; it<n; it++) {
                        matrix[i][it]=matrix[step][it];
                    }
                    for (int it=0; it<n; it++) {
                        matrix[step][it]=temp[it];
                    }
                    // pivoting done. now gaussian kora need
                    free(temp);
                    break;
                }
            }
        }

        // gaussian elimination here
        for (int i=step+1; i<n; i++) {
            double factor= matrix[i][step]/matrix[step][step];
            for(int j=step; j<n; j++) {
                if (j==step) {
                    matrix[i][j]= factor;
                }
                else {
                    matrix[i][j]+= -factor*matrix[step][j];
                }
            }
        }
    }

    if (Abs(matrix[n-1][n-1])<1e-10) {
        return 0.0;
    }

    double det=1;
    for (int i=0; i<n; i++) {
        det*=matrix[i][i];
    }

    return det*pow(-1,swap);
}

static double max_row_sum(double **matrix, const int n) {
    double max=-10;
    for (int i=0; i<n; i++) {
        double s=0;
        for (int j=0; j<n; j++) {
            s+=Abs(matrix[i][j]);
        }
        if (s>max)
            max=s;
    }
    return max;
}

static double condition_num(double **matrix, const int n) {
    double **soln= malloc(sizeof(double*)*n);
    if (!soln) {
        return NAN;
    }
    for (int i=0; i<n; i++) {
        soln[i]=malloc(sizeof(double)*n);
        if (!soln[i]) {
            free(soln);
            return NAN;
        }
    }

    inverse_matrix_without_file(matrix,soln,n);

    double val1= max_row_sum(matrix,n);
    double val2= max_row_sum(soln,n);

    // free
    for (int i=0; i<n; i++) {
        free(soln[i]);
    }
    free(soln);
    return val1*val2;
}


int main() {
    // FILE *fp= fopen("C:/Users/Rifat/Downloads/5x5_matrix.txt","r");
    // if (!fp) {
    //     return 1;
    // }

    // FILE *fo= fopen("C:/Users/Rifat/Downloads/50x50_inv_matrix.txt","w");
    // if (!fo) {
    //     fclose(fp);
    //     return 1;
    // }
    //
    // int n= total_row(fp);
    // inverse_matrix(fp,fo);
    // printf("%.6lf\n",determinant(fp));

    double mat[5][5]={
        {8.20, -1.50,0, 2.10, -0.80},
        {-2.40, 9.50, 3.10, 0, -1.90},
        {0, 1.80, 11.40, -3.20, 1.25},
        {1.75, 0, -2.60, 8.90, -1.40},
        {-1.10, 2.30, -1.50, 0, 10.60}
    };

    double **matrix= malloc(sizeof(double*)*5);
    for (int i=0; i<5; i++) {
        matrix[i]=malloc(sizeof(double)*5);
    }

    for (int i=0; i<5; i++) {
        for (int j=0; j<5; j++) {
            matrix[i][j]=mat[i][j];
        }
    }
    // printf("%.6lf\n",determinant_without_file(matrix,5));
    printf("%.10lf\n",condition_num(matrix,5));

    // fclose(fp);
    // fclose(fo);
    for (int i=0; i<5; i++) {
        free(matrix[i]);
    }
    free(matrix);
}