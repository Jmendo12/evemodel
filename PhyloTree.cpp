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
	node *parent;
	node *left;
	node *right;
	double branchLength;
	int regime;
	int numIndividuals;
	int level;
};

class PhyloTree
{
	public:
		PhyloTree(int numGenes, int verbose, std::string *treeFile, std::string *individualsFile);
		~PhyloTree();

		void insert(double key);
		node *search(double key);
		node *postOrderTraverse();
		void destroyTree();

		int readTree(std::string *fileName);
		int readNIndividuals(std::string *fileName);

	private:
		int totalIndividuals;
		int numGenes;
		int numSpecies;
		int verbose;
		int nodeStack;

		std::string *treeFile;
		std::string *individualsFile;
		
		std::vector <std::vector <double> > expressionData;
		std::vector <std::vector <double> > technicalVariance;
		
		node *root;

		void destroyTree(node *leaf);
		void insert(double key, node *leaf);
		node *search(double key, node *leaf);
};
	
// Constructor for our trees
PhyloTree::PhyloTree(int numGenes,int verbose, std::string *treeFile, std::string *individualsFile)
{
	root = NULL;

	// Assign the values to our member variables
	this->verbose = verbose;
	this->treeFile = treeFile;
	this->individualsFile = individualsFile;
	this->numSpecies = this->readTree(this->treeFile);
	this->totalIndividuals = this->readNIndividuals(this->individualsFile);
	
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


/*A basic insertion function for our binary tree */
void PhyloTree::insert(double key, node *leaf)
{
	if(key < leaf->branchLength)
	{
		if(leaf->left != NULL)
			insert(key, leaf->left);
		else
		{
			leaf->left = new node;
			leaf->left->branchLength = key;
			/* When creating a new node create a left and right child for it, but set them to null
			 * because they currently hold nothing and/or may not be needed. Also set the parent
			 * for the node as the leaf from which it came
			 */
			leaf->left->parent = leaf;
			leaf->left->left = NULL;
			leaf->left->right = NULL;
		}
	}
	else if(key >= leaf->branchLength)
	{
		if(leaf->right != NULL)
			insert(key, leaf->right);
		else
		{
			leaf->right = new node;
			leaf->right->branchLength = key;
			leaf->right->parent = leaf;
			leaf->right->left = NULL;
			leaf->right->right = NULL;
		}
	}
}

/* A function for recursively searching for a specific node from our tree */
node *PhyloTree::search(double key, node *leaf)
{
	if(leaf != NULL)
	{
		if (key == leaf->branchLength)
			return leaf;
		if(key < leaf->branchLength)
			return search(key, leaf->left);
		else
			return search(key, leaf->right);
	}
	else return NULL;
}

/* This is the public version of the insert function. It handles the case in which the root is not yet initialized
 * by the allocating memory for it and setting both child nodes to null. If the root exists
 * the recursive insert function is called
 */
void PhyloTree::insert(double key)
{
	if(root != NULL)
		insert(key, root);
	else
	{
		root = new node;
		root->branchLength = key;
		root->left = NULL;
		root->right = NULL;
	}
}

/* This is the public version of the search function. It begins with the root node */
node *PhyloTree::search(double key)
{
	return search(key, root);
}

/* The following function is a post order traversal of our trees. It will return nodes starting from
 * the left subtree, right subtree, and then the root. For example, consider the following tree
 * 5    3   1
 *  \    \  /
 *   \    2
 *    \  /
 *     4
 * The order of traversal will be 1,3,2,5,4
 */
node *PhyloTree::postOrderTraverse()
{
	postOrderTraverse(this->root);
}

/*This is the public version of the destroy tree function */
void PhyloTree::destroyTree()
{
	destroyTree(root);
}

int PhyloTree::readTree(std::string *fileName)
{
	/* Open our file represented by the fileName variable. Note that
	ifstream's constructor works only with c strings
	*/
	std::ifstream treeFile(fileName->c_str());
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
		should be a new line before the EOF, or else throw an exception,
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

		/*for(int i = 0; i < numSpecies; i++)
		{
			this->root[i].up[0] = -1;
			this->root[i].up[1] = -1;
		}

		int root = 2 * numSpecies - 2;
		this->nodeStack = root;
		this->root[root].down = -1;
		*/
		//getclade();

		treeFile.close();

		//addLevelsRecurse(root, 0);
		
	}

	return numSpecies;
}

int PhyloTree::readNIndividuals(std::string *fileName)
{
	int totalIndividuals = 0;
	std::ifstream treeFile(fileName->c_str());

	if(this->verbose > 1)
		std::cout << "Reading number of individuals file" << std::endl;

	if(!treeFile.is_open())
	{
		std::cout << "Cannon open number of individuals file" << std::endl;
		throw std::exception();
	}

	/*for(int i = 0; i < this->numSpecies; i++)
	{
		treefile >> this->node[i]->numIndividuals;
		totalIndividuals += this->node[i]->numIndividuals;
	}*/
	
	treeFile.close();
	return totalIndividuals;	
}
