#include <iostream>
#include "omp.h"
#include <ctime>
#include <cstdlib>
#include <string>
#include <iomanip>
#include <chrono>
#include <tuple>
#include "math.h"
#include <fstream>
using namespace std;

int countwords(string mytemplate)
{
	ifstream fin;
	fin.open(mytemplate);
	string word;
	int count=0;
	while(!fin.eof())
	{
		fin>>word;
		count++;
	}
	fin.close();
	return count;
}

int main() {
//	*** INIT ***

	string line;
	string textFull;
	string mytemplate("MALCOLM");

	string filename = "mytext.txt";
	cout << "Number of words " << countwords(filename) << endl;

	ifstream mytext(filename);

        if (mytext.is_open()) {
		while (getline(mytext, line)) {
			textFull += line;
			textFull += " ";
 		}
                mytext.close();

        } else {
                cout << "File isn't open";
	}

	unsigned int threadsN = 8;
	omp_set_num_threads(threadsN);

	ifstream wordsStream(filename);


//	*** FIND ***

//	*** PARALLEL ***
	unsigned int sectionN = threadsN;
	unsigned int lengthSection = textFull.length() / sectionN;
	unsigned int lengthTail = textFull.length() % sectionN;

	auto begin = chrono::steady_clock::now();
	unsigned int s = 0;
	# pragma omp parallel for reduction(+:s)
	for (int sectionNumber = 0; sectionNumber < sectionN; sectionNumber++) {
		unsigned int localPos = sectionNumber*lengthSection;
		string textLocal(textFull.substr(localPos, lengthSection + mytemplate.length()-1));
		int foundPos = -1;
		bool first = true;
		while (foundPos!=string::npos || first) {
			first = false;
			foundPos = textLocal.find(mytemplate, foundPos + 1);
			if (foundPos!=string::npos)
				s+=1;
		}
	}

	unsigned int localPos = sectionN*lengthSection;
	string textLocal(textFull.substr(localPos, lengthTail + mytemplate.length()-1));
	int foundPos = -1;
	bool first = true;

	while (foundPos!=string::npos || first) {
		first = false;
		foundPos = textLocal.find(mytemplate, foundPos + 1);
		if (foundPos!=string::npos) {
			s += 1;
		}
	}
        auto end = chrono::steady_clock::now();
        auto t1 =
		chrono::duration_cast<chrono::microseconds>(end - begin).count();

//	*** SERIAL ***
	begin = chrono::steady_clock::now();
	unsigned int s2 = 0;
	foundPos = -1;
	first = true;
	while (foundPos!=string::npos || first) {
		first = false;
		foundPos = textFull.find(mytemplate, foundPos + 1);
		if (foundPos!=string::npos) {
			s2 += 1;
		}
	}
        end = chrono::steady_clock::now();
	auto t2 =
		chrono::duration_cast<chrono::microseconds>(end - begin).count();

	cout << "Parallel " <<  t1 << ", Serial " << t2 << endl;
}
