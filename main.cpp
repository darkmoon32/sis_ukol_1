#include <math.h>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

double evalFunction(double num)
{
    return ((pow(num, 3)) * (sqrt(1 + (num / 2))));
}

double obdelnikove(double h, double dolniMez, double horniMez)
{
    double vysledek = 0;
    while(dolniMez < horniMez)
    {
        vysledek += h * evalFunction(dolniMez + h / 2);
        dolniMez += h;
    }
    return vysledek;
}

double lichobeznikove(double h, double dolniMez, double horniMez)
{
    double vysledek = 0;

    while(dolniMez < horniMez)
    {
        vysledek += (h / 2) * (evalFunction(dolniMez) + evalFunction(dolniMez + h));
        dolniMez += h;
    }
    return vysledek;
}

double simpsonova(double h, double dolniMez, double horniMez)
{
    double vysledek = 0;

    while(dolniMez < horniMez)
    {
        vysledek += (h / 6) * (evalFunction(dolniMez) + 4 * evalFunction(dolniMez + h) + evalFunction(dolniMez + 2 * h));
        dolniMez += h;
    }
    return vysledek;
}

std::vector<double> * romberg(double a, double b, double h)
{
    std::vector<double> *results = new std::vector<double>(6);

    results->at(0) = obdelnikove(h, a, b);      // t00
    results->at(1) = obdelnikove(h / 2, a, b);  // t01
    results->at(3) = obdelnikove(h / 4, a, b);  // t02

    results->at(2) = (1 / (pow(4, 1) - 1)) * (pow(4, 1) * results->at(1) - results->at(0)); // t10
    results->at(4) = (1 / (pow(4, 1) - 1)) * (pow(4, 1) * results->at(3) - results->at(1)); // t11
    results->at(5) = (1 / (pow(4, 2) - 1)) * (pow(4, 2) * results->at(4) - results->at(2)); // t20

    return results;
}

double firstDerivTwoPoint(double h, double arg){
    return (1 / h) * (evalFunction(arg) - evalFunction(arg - h));
}

double firstDerivThreePoint(double h, double arg){
    return (1 / (2 * h)) * (-3 * evalFunction(arg) + 4 * evalFunction(arg + h) - evalFunction(arg + 2 * h));
}

double firstDerivFourPoint(double h, double arg){
    return (1 / (6 * h)) * (-11 * evalFunction(arg) + 18 * evalFunction(arg + h) - 9 * evalFunction(arg + 2 * h) + 2 * evalFunction(arg + 3 * h));
}

double firstDerivFivePoint(double h, double arg){
    return (1 / (12 * h)) * (-25 * evalFunction(arg) + 48 * evalFunction(arg + h) - 36 * evalFunction(arg + 2 * h) + 16 * evalFunction(arg + 3 * h) - 3 * evalFunction(arg + 4 * h));
}

double secondDerivThreePoint(double h, double arg){
    return (1 / pow(h, 2)) * (evalFunction(arg) - 2 * evalFunction(arg + h) + evalFunction(arg + 2 * h));
}

double secondDerivFourPoint(double h, double arg){
    return (1 / pow(h, 2)) * (2 * evalFunction(arg) - 5 * evalFunction(arg + h) + 4 * evalFunction(arg + 2 * h) - evalFunction(arg + 3 * h));
}

std::vector<double> * richardson(double h, double arg, double (*pFunc)(double, double))
{
    std::vector<double> *results = new std::vector<double>(6);

    results->at(0) = pFunc(h, arg);      // d00
    results->at(1) = pFunc(h / 2, arg);  // d01
    results->at(3) = pFunc(h / 4, arg);  // d02

    results->at(2) = (1 / (pow(4, 1) - 1)) * (pow(4, 1) * results->at(1) - results->at(0)); // d10
    results->at(4) = (1 / (pow(4, 1) - 1)) * (pow(4, 1) * results->at(3) - results->at(1)); // d11
    results->at(5) = (1 / (pow(4, 2) - 1)) * (pow(4, 2) * results->at(4) - results->at(2)); // d20

    return results;
}

int main(int argc, char *argv[])
{
    ofstream myfile("output.csv");
    if(myfile){
        myfile.close();
        remove("output.csv");
    }
    myfile.open("output.csv");

    //3.73532
    std::vector<double> *rombergResult;
    myfile << "Interval;Obdelnik;Lichobeznik;Simpson;Analaticky\n";
    for(int i = 0 ; i < 20 ; i++)
    {
        myfile << std::to_string(4 / pow(2 , i)) + ';' + to_string(obdelnikove(4 / pow(2, i), -2, 2)) + ';' +
                  to_string(lichobeznikove(4 / pow(2, i), -2, 2)) + ';' + to_string(simpsonova(4 / pow(2, i), -2, 2)) + ';' + to_string(3.73532) + '\n';
    }

    rombergResult = romberg(-2, 2, 2);
    myfile << "\nT;Zpresneni k=1;Zpresneni k=2\n";
    myfile << to_string(rombergResult->at(0)) + '\n';
    myfile << to_string(rombergResult->at(1)) + ';' + to_string(rombergResult->at(2)) + '\n';
    myfile << to_string(rombergResult->at(3)) + ';' + to_string(rombergResult->at(4)) + ';' + to_string(rombergResult->at(5)) + '\n';

    myfile << "Interval;Dvou;Tri;Ctyr;Peti;Analiticky;Tri;Ctyr;Analiticky\n";
    for(int i = 0 ; i < 20 ; i++){
        myfile << to_string(4 / pow(2 , i)) + ';' + to_string(firstDerivTwoPoint(4 / pow(2, i), -1)) + ';' + to_string(firstDerivThreePoint(4 / pow(2, i), -1)) + ';' +
                  to_string(firstDerivFourPoint(4 / pow(2, i), -1)) + ';' + to_string(firstDerivFivePoint(4 / pow(2, i), -1)) + ';' + to_string(1.76777) + ';' +
                  to_string(secondDerivThreePoint(4 / pow(2, i), -1)) + ';' + to_string(secondDerivFourPoint(4 / pow(2, i), -1)) + ';' + to_string(-1.94454) + '\n';
    }

    rombergResult = richardson(0.25, -1, &firstDerivFourPoint);
    myfile << "\nT;Zpresneni k=1;Zpresneni k=2\n";
    myfile << to_string(rombergResult->at(0)) + '\n';
    myfile << to_string(rombergResult->at(1)) + ';' + to_string(rombergResult->at(2)) + '\n';
    myfile << to_string(rombergResult->at(3)) + ';' + to_string(rombergResult->at(4)) + ';' + to_string(rombergResult->at(5)) + '\n';

    rombergResult = richardson(0.25, -1, &secondDerivFourPoint);
    myfile << "\nT;Zpresneni k=1;Zpresneni k=2\n";
    myfile << to_string(rombergResult->at(0)) + '\n';
    myfile << to_string(rombergResult->at(1)) + ';' + to_string(rombergResult->at(2)) + '\n';
    myfile << to_string(rombergResult->at(3)) + ';' + to_string(rombergResult->at(4)) + ';' + to_string(rombergResult->at(5)) + '\n';
    rombergResult->~vector();

    myfile.close();
    return 0;
}
