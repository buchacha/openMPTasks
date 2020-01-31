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

int main() {
//	*** INIT ***

	string line;
	string textFull;
	string mytemplate;

	string filename = "mytext.txt";
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

//	*** TEST ***
	if (wordsStream.is_open()) {
		while (!wordsStream.eof()) {
			wordsStream >> mytemplate;
			//	*** FIND ***

			unsigned int sectionN = threadsN;
			unsigned int lengthSection = textFull.length() / sectionN;
			unsigned int lengthTail = textFull.length() % sectionN;

			unsigned int s = 0;
			# pragma omp parallel for reduction(+:s)
			for (int sectionNumber = 0; sectionNumber < sectionN; sectionNumber++) {
				unsigned int localPos = sectionNumber*lengthSection;
				string textLocal(textFull.substr(localPos, lengthSection + mytemplate.length()-1));
				int foundPos = -1;
				bool first = true;
				# pragma omp critical
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
			cout << s << " == " << s2 << " " << mytemplate << endl;
			if (s != s2) {
				cout << "sectionN " << sectionN << endl;
				cout << "lengthSection " << lengthSection << endl;
				cout << "lengthTail " << lengthTail << endl;
				cout << "s " << s << endl;
				cout << "s2 " << s2 << endl;
				cout << "temp " << mytemplate << endl;
			}
		}
		wordsStream.close();
	}
}
