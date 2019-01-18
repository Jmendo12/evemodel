#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

void printVectorDoubles(std::vector<double> *v);
void print2DVectorDoubles(std::vector< std::vector<double>> *v);

void writeVectorDoubles(std::vector<double> *v, std::string *fileName);
void writeLRTVectors(std::vector<double> *nulls, std::vector<double> *alts, std::string *fileName);
void write2DVectorDoubles(std::vector< std::vector<double>> *v, std::string *fileName);
void write2DVectorDoublesAnnotated(std::vector< std::vector<double>> *v, std::vector<std::string> *rowNames, 
		std::string *fileName);
void write3DVectorDoubles(std::vector< std::vector < std::vector<double>>> *v, std::string *fileName);
void write4DVectorDoubles(std::vector< std::vector < std::vector < std::vector<double>>>> *v, std::string *fileName);
void writeStringToFile(std::string *toWrite, std::string *fileName);

void readStringFromFile(std::string *destString, std::string *fileName, int numCharacters);
void readRowNames2DVectorDoubles(std::vector< std::vector <double>> *v, int numCols, 
		std::vector<std::string *lineNames, std::string *fileName);
void read4DVectorDoubles(std::vector <std::vector <std::vector <std::vector <double>>>> *v, std::string *fileName);
void read4DVectorDoublesto3DVectorDoubles(std::vector < std::vector < std::vector <double>>> *v, std::string *fileName);
