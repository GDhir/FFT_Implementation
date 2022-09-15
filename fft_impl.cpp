#include<iostream>
#include<vector>
#include<cmath>
#include <complex>
#include "gnuplot-iostream.h"
#define DJ_FFT_IMPLEMENTATION // define this in exactly *one* .cpp file
#include "dj_fft.h"

using namespace std;

using complext = complex<double>;

void fft_alg( vector<complext>& input, vector<complext>& output, int i, int j, int stride, int N ) {

    if( i == j ) {
        output[0] = input[i];
        return;
    }

    int halfN = N/2;

    vector<complext> oddoutput( halfN, complext(0, 0) );
    vector<complext> evenoutput( halfN, complext(0, 0) );;

    fft_alg( input, evenoutput, i, j - stride, stride*2, halfN );
    fft_alg( input, oddoutput, i + stride, j, stride*2, halfN );

    for( int idx = 0; idx < N; idx++ ) {

        output[ idx ] = evenoutput[idx%halfN] + oddoutput[idx%halfN]*complext( cos( ( 2*M_PI* idx )/N ), -sin( ( 2*M_PI* idx )/N ) );

    }

}

// void demo_basic() {
//     Gnuplot gp;
//     // For debugging or manual editing of commands:
//     //Gnuplot gp(std::fopen("plot.gnu", "w"));
//     // or
//     //Gnuplot gp("tee plot.gnu | gnuplot -persist");
//     gp << "set terminal qt\n";

//     std::vector<std::pair<double, double>> xy_pts_A;
//     for(double x=-2; x<2; x+=0.01) {
//         double y = x*x*x;
//         xy_pts_A.emplace_back(x, y);
//     }

//     std::vector<std::pair<double, double>> xy_pts_B;
//     for(double alpha=0; alpha<1; alpha+=1.0/24.0) {
//         double theta = alpha*2.0*3.14159;
//         xy_pts_B.emplace_back(cos(theta), sin(theta));
//     }

//     gp << "set xrange [-2:2]\nset yrange [-2:2]\n";
//     gp << "plot '-' with lines title 'cubic', '-' with points title 'circle'\n";
//     gp.send1d(xy_pts_A);
//     gp.send1d(xy_pts_B);

//     //pause_if_needed();
// }

void write_output( const string& filename, const vector<complext>& input, const vector<complext>& output ) {

     std::ofstream myFile(filename);
    
    // Send data to the stream
    for(size_t i = 0; i < input.size(); ++i)
    {
        myFile << input.at(i).real() << "," << input.at(i).imag() << "," << output.at(i).real() << "," << output.at(i).imag() << "\n";
    }
    
    myFile << "Magnitude" <<"\n";

    for(size_t i = 0; i < input.size(); ++i)
    {
        myFile << sqrt( pow( output.at(i).real(), 2 ) + pow( output.at(i).imag(), 2 ) ) << "\n";
    }

    // Close the file
    myFile.close();



}

int main() {

    int N = 32;
    vector<complext> input(N, complext(0, 0));
    vector<complext> output(N, complext(0, 0));

    //input[15] = make_pair(1, 0);
    for( int i = 0; i < N; i++ ) {

        input[i] = complext( sin( ( 2*M_PI/N )*i ), 0 );

    }

    fft_alg( input, output, 0, N - 1, 1, N );

    for( auto& elem: output ) {
        elem /= sqrt(N);
    }

    auto fftExact = dj::fft1d(input, dj::fft_dir::DIR_FWD);

    write_output("currentFFTFile.csv", input, output);

    write_output( "exactFFTFile.csv",input, fftExact );

    return 0;

}
