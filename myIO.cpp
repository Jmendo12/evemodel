#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

void printVectorDoubles(std::vector<double> *v)
{
	for(int i = 0; i < v->size(); i++)
		std::cout << v->at(i) << " ";
	std::cout << std::endl;
}

void print2DVectorDoubles(std::vector< std::vector <double> > *v)
{
	for(int i = 0; i < v->size(); i++)
	{
		for(int j = 0; j < v->at(i).size(); j++)
			std::cout << v->at(i).at(j) << " ";

		std::cout << std::endl;
	}
	std::cout << std::endl;
}

void writeVectorDoubles(std::vector<double> *v, std::string *fileName)
{
	std::ofstream writeFile (fileName->c_str());

	if(writeFile.is_open())
	{
		for(int i = 0; i < v->size(); i++)
			writeFile << v->at(i) << " ";
		
		std::cout << std::endl;
		writeFile.close();
	}
	else
		std::cout << "Unable to open the file for printing your vector of doubles." << std::endl;
}

void writeLRTVectors(std::vector<double> *nulls, std::vector<double> *alts, std::string *fileName)
{
	std::ofstream writeFile (fileName->c_str());

	if(writeFile.is_open())
	{
		for(int i = 0; i < nulls->size() || i < alts->size(); i++)
			writeFile << 2.0 * alts->at(i) - 2.0 * nulls->at(i) << " ";

		std::cout << std::endl;

		writeFile.close();
	}
	else
		std::cout << "Unable to open the file for printing your LRT vectors." << std::endl;
}

void write2DVectorDoubles(std::vector< std::vector <double> > *v, std::string *fileName)
{
	std::ofstream writeFile (fileName->c_str());

	if(writeFile.is_open())
	{
		for(int i = 0; i < v->size(); i++)
		{
			for(int j = 0; j < v->at(i).size(); j++)
				writeFile << v->at(i).at(j) << " ";

			std::cout << std::endl;
		}

		writeFile.close();
	}
	else
		std::cout << "Unable to open the file for printing your 2D vector of doubles." << std::endl;
}

void write2DVectorDoublesAnnotated(std::vector< std::vector <double> > *v, std::vector<std::string> *rowNames, 
		std::string *fileName)
{
	std::ofstream writeFile (fileName->c_str());

	if(writeFile.is_open())
	{
		if(v->size() == rowNames->size())
		{
			std::cout << "Number of rows: " << v->size()  << std::endl;
			
			for(int i = 0; i < v->size(); i++)
			{
				std::cout << "Row name: " << rowNames->at(i);
				for(int j = 0; j < v->at(i).size(); j++)
					writeFile << v->at(i).at(j) << " ";

				std::cout << std::endl;
			}
		}
		else
			std::cout << "The amount of names for rows differs from the amount of rows, they should be"
				<< " the same." << std::endl;
		writeFile.close();
	}
	else
		std::cout << "Unable to open the file printing your annotated 2D vector of doubles" << std::endl;
}

void write3DVectorDoubles(std::vector < std::vector < std::vector <double> > > *v, std::string *fileName)
{
	std::ofstream writeFile (fileName->c_str());

	if(writeFile.is_open())
	{
		for(int i = 0; i < v->size(); i++)
		{
			for(int j = 0; j < v->at(i).size(); j++)
			{
				for(int k = 0; k < v->at(i).at(j).size(); k++)
					writeFile << v->at(i).at(j).at(k) << " ";

				std::cout << std::endl;
			}
			std::cout << std::endl;
		}

		writeFile.close();
	}
	else
		std::cout << "Unable to open the file for printing your 3D double vector." << std::endl;
}

void write4DVectorDoubles(std::vector < std::vector < std::vector < std::vector <double> >>> *v, std::string *fileName)
{
	std::ofstream writeFile(fileName->c_str());

	if(writeFile.is_open())
	{
		for(int i = 0; i < v->size(); i++)
		{
			for(int j = 0; j < v->at(i).size(); j++)
			{
				for(int k = 0; k < v->at(i).at(j).size(); k++)
				{
					for(int l = 0; l < v->at(i).at(j).at(k).size(); l++)
						writeFile << v->at(i).at(j).at(k).at(l) << " ";

					std::cout << std::endl;
				}
				std::cout << std::endl;
			}
			std::cout << std::endl;
		}
		writeFile.close();
	}
	else
		std::cout << "Unable to open the file for printing your 4D double vector." << std::endl;
}

void writeStringToFile(std::string *toWrite, std::string *fileName)
{
	std::ofstream writeFile(fileName->c_str());

	if(writeFile.is_open())
	{
		writeFile << *toWrite;
		writeFile.close();
	}
	else
		std::cout << "Unable to open the file for printing a string." << std::endl;
}

void readStringFromFile(std::string *destString, std::string *fileName, int numCharacters)
{
	std::ifstream readFile(fileName->c_str());

	if(readFile.is_open())
	{
		for(int i = 0; i < numCharacters && !readFile.eof(); i++)
		{
			destString->push_back(readFile.get());
		}

		readFile.close();
	}
	else
		std::cout << "Unable to open the file for reading a string. " << std::endl;
}

void readRowNames2DVectorDoubles(std::vector < std::vector <double> > *v, int numCols, 
		std::vector<std::string> *lineNames, std::string *fileName)
{
	std::ifstream readFile(fileName->c_str());
	std::string rowStream;
	int numRows = 0;

	if(readFile.is_open())
	{
		// Get the first line of the file and place it into a string
		getline(readFile, rowStream);
		//Convert the string to an integer by creating a stringstream and streaming it into numRows	
		std::stringstream geek(rowStream);
		geek >> numRows;
			
		for(int i = 0; i < numRows; i++)
		{
			/* Get the line until a space is found and place the result in rowStream.
			 * Then add rowStream to the back of the vector
			 */
			getline(readFile, rowStream, ' ');
			lineNames->push_back(rowStream);

			for(int j = 0; j < numCols; j++)
				readFile >> v->at(i)[j];
		}

		readFile.close();
	}
	else
		std::cout << "Unable to open the file for reading a 2D vector of doubles with row names." << std::endl;
}

//Assume that we know how the dimension of the vector when we create the vector
void read4DVectorDoubles(std::vector < std::vector < std::vector < std::vector <double>>>> *v, std::string *fileName)
{
	std::ifstream readFile(fileName->c_str());
	if(readFile.is_open())
	{
		for(int i = 0; i < v->size(); i++)
		{
			for(int j = 0; j < v->at(i).size(); j++)
			{
				for(int k = 0; k < v->at(i).at(j).size(); k++)
				{
					for(int l = 0; l < v->at(i).at(j).at(k).size(); l++)
					{
						readFile >> v->at(i).at(j).at(k)[l];
					}
				}
			}
		}
		readFile.close();
	}
	else
		std::cout << "Unable to open the file for reading a 4D vector of doubles." << std::endl;
}

//Assume we know the dimension of the vector we pass into as an argument when we create the vector
void read4DVectorDoublesTo3DVectorDoubles(std::vector < std::vector < std::vector <double>>> *v, std::string *fileName)
{
	std::ifstream readFile(fileName->c_str());
	if(readFile.is_open())
	{
		for(int i = 0; i < v->size(); i++)
		{
			for(int j = 0; j < v->at(i).size(); j++)
			{
				for(int k = 0; k < v->at(i).at(j).size(); k++)
					readFile >> v->at(i).at(j)[k];
			}
		}
		readFile.close();
	}
	else
		std::cout << "Unable to open the file for reading a 4D vector of doubles." << std::endl;
}


