#include <iostream>
#include <fstream>
#include <cmath>
#include <cctype>
#include <ctime>
#include <vector>
#include <string>
#include <exception>

//#include "myio.h" not yet implemented

struct node
{
	int down;
	int up[2];
	double bl;
	node *left;
	node *right;
	int regime;
	int numIndividuals;
	int level;
};

class PhyloTree
{
	public:
		PhyloTree(int numSpecies, int numGenes, int totalIndividuals, int verbose);
		~PhyloTree();

		void insert(int key);
		node *search(int key);
		void destroyTree();

		int readTree(std::string fileName);


	private:
		int totalIndividuals;
		int numGenes;
		int numSpecies;
		int verbose;
		int nodeStack;
		
		std::vector <std::vector <double> > expressionData;
		std::vector <std::vector <double> > technicalVariance;
		
		node *root;

		void destroyTree(node *leaf);
		void insert(int key, node *leaf);
		node *search(int key, node *leaf);
};
	
	// Constructor for our trees
PhyloTree::PhyloTree(int numSpecies, int numGenes, int totalIndividuals, int verbose)
{
	root = NULL;

	// Assign the values to numSpecies, numGenes, and verbose
	this->numSpecies = numSpecies;
	this->numGenes = numGenes;
	this->totalIndividuals = totalIndividuals;
	this->verbose = verbose;
	/* expressionData is initialized as a vector with numGenes vectors as rows and 
	totalIndividuals doubles as columns within those rows */
	std::vector <std::vector <double> > expressionData (numGenes, std::vector <double> (totalIndividuals)); 

	/* technical Variance is initialized as a vector with numGenes vectors as rows and
	totalIndividuals doubles as columns within those rows. The values for all indices are then set to 0
	*/
	std::vector <std::vector <double> > technicalVariance (numGenes, std::vector <double> (totalIndividuals, 0.0));
}

/* The destructor for our trees; since our vectors are not dynamically allocated,
we need not destroy them */
PhyloTree::~PhyloTree()
{
	destroyTree();
}

/* The recursive destructor that does the work involved with deletion. If the
leaf is not null then we delete its possible children, and delete the leaf */
void PhyloTree::destroyTree(node *leaf)
{
	if(leaf != NULL)
	{
		destroyTree(leaf->left);
		destroyTree(leaf->right);
		delete leaf;
	}
}

int PhyloTree::readTree(std::string fileName)
{
	/* Open our file represented by the fileName variable. Note that
	ifstream's constructor works only with c strings
	*/
	std::ifstream treeFile(fileName.c_str());
	int numSpecies;
	char c;

	if(this->verbose > 1)
		std::cout << "Reading tree (with species names as integers) " << std::endl;

	/* If the file does not exist or cannot be opened throw an exception;
	note that this exception should be caught in main, and the program should
	be subsequently terminated via EXIT_FAILURE */
	if(!treeFile.is_open())
	{
		std::cout << "Cannot open tree file: " << fileName << "!" << std::endl;
		throw std::exception(); 
	}
	else
	{
		treeFile >> numSpecies;

		if(this->verbose > 1)
			std::cout << "Number of species: " << numSpecies << std::endl;

		if(numSpecies < 4)
		{
			std::cout << "At least four species are needed for this program" << std::endl;
			throw std::exception();
		}

		/*Loop to ensure that our input file is of the format we need. There
		should be a new line before the EOF, or else thrown an exception,
		which should be caught in main */
		while(treeFile.get(c))
		{
			if(treeFile.eof())
			{
				std::cout << "Error reading treeFile" << std::endl;
				throw std::exception();
			}
			else if(c == '\n')
				break;
		}

		for(int i = 0; i < numSpecies; i++)
		{
			this->root[i].up[0] = -1;
			this->root[i].up[1] = -1;
		}

		int root = 2 * numSpecies - 2;
		this->nodeStack = root;
		this->root[root].down = -1;
		//getclade();

		treeFile.close();

		//addLevelsRecurse(root, 0);
	}

	return numSpecies;
}