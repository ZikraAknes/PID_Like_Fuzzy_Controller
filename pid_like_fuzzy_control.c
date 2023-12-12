#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <dir.h>
#include "libraries/pbPlots.c"
#include "libraries/supportLib.c"

#define abs_float(x) ((x)>0?(x):-(x))

// Global Variables
char rules[3][2];
char rules_list[3] = {'N', 'Z', 'P'};
double m[3][2];
double w[8];
double force[8];

// Function to add char in a string
char add(char *s, const int len, char x, int pos){
    s[len] = 0;
    memmove(s+pos+1, s+pos, len-pos+1);
    s[pos]=x;
}

// Function to find the index from rules_list
int find_idx(char x){
    int idx;
    for(idx = 0; idx<5; idx++){
        if(rules_list[idx] == x){
            break;
        }
    }

    return idx;
}

// Calculations for rise time
double rise_time(double time[], double y[], int iteration){
    double y_target = 0.1;
    double ten_percent;
    double ninety_percent;

    for(int i = 0; i<2; i++){
        int idx;
        for(idx=0; idx<iteration; idx++){
            if(y[idx] > y_target){
                break;
            }
        }

        double xtr1 = time[idx-1];
        double xtr2 = time[idx];
        double ytr1 = y[idx-1];
        double ytr2 = y[idx];
        double grad = (ytr2 - ytr1) / (xtr2 - xtr1);
        double y_intercept = ytr1-(grad*xtr1);

        if(i == 0)
            ten_percent = ((y_target-y_intercept)/grad);
        else
            ninety_percent = ((y_target-y_intercept)/grad);

        y_target = 0.9;
    }

    return (ninety_percent - ten_percent);
}

// Calculations for steady state error
double steady_state_error(double y[], double d[], int iteration){
    double sse = INFINITY;
    for(int i = 1; i<iteration; i++){
        if(abs_float(y[i] - y[i-1]) < 1e-5){
            sse = d[i] - y[i];
            break;
        }
    }

    return sse;
}

// Calculations for overshoot
double overshoot(double y[], double d[], int iteration){
    double os = 0;
    for(int i = 0; i<iteration; i++){
        if(y[i] > os && y[i] > d[i]){
            os = y[i];
        }
    }

    if(os != 0)
        os -= d[0];

    return os;
}

// Function to calculate mu
void calculate_mu(double error, double a, double b, int idx){
    double m1 = (b - error)/(b-a);
    double m2 = (error - a)/(b-a);

    m[idx][0] = m1;
    m[idx][1] = m2;
}

// Find rules and calculate mu
void membership_function(double error, double e[], int idx, double h){
    double e1 = e[0];
    double e2 = e[1];
    double e3 = e[2];

    if(h != 1){
        e2 = h * error;
    }
    
    if(error <= e1){
        m[idx][0] = 1;
        m[idx][1] = 0;
        rules[idx][0] = 'N';
        rules[idx][1] = 'Z';
    }else if(error < e2){
        calculate_mu(error, e1, e2, idx);
        rules[idx][0] = 'N';
        rules[idx][1] = 'Z';
    }else if(error == e2){
        m[idx][0] = 0;
        m[idx][1] = 1;
        rules[idx][0] = 'N';
        rules[idx][1] = 'Z';
    }else if(error < e3){
        calculate_mu(error, e2, e3, idx);
        rules[idx][0] = 'Z';
        rules[idx][1] = 'P';
    }else if(error >= e3){
        m[idx][0] = 0;
        m[idx][1] = 1;
        rules[idx][0] = 'Z';
        rules[idx][1] = 'P';
    }
}

void fuzzy_control(double defuzz[], double err_range[], double delta_err_range[], double double_delta_range[], int num_plant, int iteration, double h){
    char path[100] = "outputs/plant.png";
    char plant[2];
    char plant_str[10] = "PLANT ";
    wchar_t plant_str2[10];

    sprintf(plant, "%d", num_plant);

    add(path, 100, *plant, 13);
    add(plant_str, 10, *plant, 6);

    mbstowcs(plant_str2, plant_str, 10);

    printf("\n=======%s=======\n", plant_str);
    
    // Output Variable
    double y[500] = {0.0};
    double x2[500] = {0.0};
    double x1[500] = {0.0};
    
    // Fuzzy output variable
    double u[500] = {0.0};
    // Desired Output
    double d[500] = {0.0};
    // Time variable
    double time[500];
    // Previous error
    double err_prev;
    double delta_err_prev;
    // Current error
    double err = d[0] - y[0];
    double delta_err;
    double double_delta_err;

    // Get all defuzzifier value
    double N = defuzz[0];
    double Z = defuzz[1];
    double P = defuzz[2];

    // Force of rules matrix
    double rules_matrix[3][3][3] = {{{N, N, N}, {N, N, N}, {P, P, P}},
                                    {{N, N, N}, {Z, Z, Z}, {P, P, P}},
                                    {{N, N, N}, {P, P, P}, {P, P, P}}};

    // Set value for time and desired output
    for(int i = 0; i<iteration; i++){
        time[i] = (double)i;
        d[i] = 1;
    }

    for(int k = 0; k<iteration-1; k++){
        // Errors calculations
        err_prev = err;
        err = d[k] - y[k];
        delta_err_prev = delta_err;
        delta_err = err - err_prev;
        double_delta_err = delta_err - delta_err_prev;

        // Fuzzifier
        membership_function(err, err_range, 0, h);
        membership_function(delta_err, delta_err_range, 1, 1);
        membership_function(double_delta_err, double_delta_range, 2, 1);

        // Calculate weights
        for(int i = 0; i<2; i++){
            for(int j = 0; j<2; j++){
                for(int z = 0; z<2; z++){
                    w[(((i*2)+j)*2) + z] = m[0][i]*m[1][j]*m[2][z];
                }
            }
        }

        // Defuzzifier
        for(int i = 0; i<2; i++){
            int i_idx = find_idx(rules[0][i]); 
            
            for(int j = 0; j<2; j++){
                int j_idx = find_idx(rules[1][j]);

                for(int z = 0; z<2; z++){
                    int z_idx = find_idx(rules[2][z]);

                    force[(((i*2)+j)*2) + z] = rules_matrix[j_idx][i_idx][z_idx];
                }
            }
        }

        u[k+1] = 0;
        double w_tot = 0;
        
        // Calculate fuzzy output
        for(int i = 0; i<8; i++){
            w_tot += w[i];
            u[k+1] += w[i]*force[i];
        }
        u[k+1] /= w_tot;

        // Fuzzy output as Δu(k)
        u[k+1] = u[k] + u[k+1];

        if(num_plant == 1){
            // 1st Plant
            y[k+1] = ((1.01*y[k]) + (0.01*(pow(y[k], 2))) + (0.03*u[k+1]));
        }else if(num_plant == 2){
            // 2nd Plant
            x1[k+1] = (x1[k]+(0.01*x2[k])+(0.01*u[k+1]));
            x2[k+1] = ((0.1*x1[k])+(0.97*x2[k]));
            y[k+1] = x1[k+1];
        }
    }

    // Calculate Rise Time
    double tr = rise_time(time, y, iteration);

    // Calculate SSE
    double sse = steady_state_error(y, d, iteration);

    // Calculate Overshoot
    double os = overshoot(y, d, iteration);

    printf("\nRise Time : %f", tr);
    printf("\nOvershoot : %f", os);
    if(sse != INFINITY)
        printf("\nSSE       : %f\n", sse);
    else
        printf("\nSSE       : %s", "System is Unsteady\n");

    // Create directory for saving plot results
    mkdir("outputs");

    // Plot output
    ScatterPlotSeries *series = GetDefaultScatterPlotSeriesSettings();
	series->xs = time;
	series->xsLength = iteration;
	series->ys = y;
	series->ysLength = iteration;
	series->linearInterpolation = true;
	series->lineType = L"solid";
	series->lineTypeLength = wcslen(series->lineType);
	series->lineThickness = 2;
	series->color = CreateRGBColor(0, 0, 1);

    // Plot desired output
    ScatterPlotSeries *series2 = GetDefaultScatterPlotSeriesSettings();
	series2->xs = time;
	series2->xsLength = iteration;
	series2->ys = d;
	series2->ysLength = iteration;
	series2->linearInterpolation = true;
	series2->lineType = L"dashed";
	series2->lineTypeLength = wcslen(series2->lineType);
	series2->lineThickness = 2;
	series2->color = GetGray(0.6);

    // Plot settings
	ScatterPlotSettings *settings = GetDefaultScatterPlotSettings();
	settings->width = 1000;
	settings->height = 600;
	settings->autoBoundaries = false;
    settings->yMax = 1.2;
    settings->xMax = iteration;
	settings->autoPadding = true;
    settings->title = plant_str2;
    settings->titleLength = wcslen(settings->title);
    settings->xLabel = L"Y(k)";
    settings->xLabelLength = wcslen(settings->xLabel);
    settings->yLabel = L"Steps";
    settings->yLabelLength = wcslen(settings->yLabel);
	ScatterPlotSeries *s [] = {series2, series};
	settings->scatterPlotSeries = s;
	settings->scatterPlotSeriesLength = 2;

	RGBABitmapImageReference *canvasReference = CreateRGBABitmapImageReference();
	DrawScatterPlotFromSettings(canvasReference, settings);

	size_t length;
	double *pngdata = ConvertToPNG(&length, canvasReference->image);
	WriteToFile(pngdata, length, path);
	DeleteImage(canvasReference->image);

    printf("\n=====================\n", plant_str);
    printf("\n[INFO] Plot Saved To: '%s'\n", path);
}

int main(){
    // Defuzzifier value
    double uk[3] = {-10, 0, 10};

    // Error initial parameter value
    double ek[3] = {-1, 0, 1};

    // Delta error initial parameter value
    double delta_ek[3] = {-1, 0, 1};

    // Double delta error initial parameter value
    double double_delta_ek[3] = {-1, 0, 1};

    // H initial value
    double h = 1;

    // Calculate and Plot 1st Plant
    delta_ek[0] = -0.5;
    delta_ek[1] = -0;
    delta_ek[2] = 0.5;
    h = 0.8;
    fuzzy_control(uk, ek, delta_ek, double_delta_ek, 1, 100, h);

    // Calculate and Plot 2nd Plant
    delta_ek[0] = -0.15;
    delta_ek[1] = -0;
    delta_ek[2] = 0.15;
    h = 0.45;
    fuzzy_control(uk, ek, delta_ek, double_delta_ek, 2, 200, h);

    printf("\n");

    return 0;
}