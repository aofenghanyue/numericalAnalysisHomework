#include "homework3.h"
#include <ctime>

int main() {
	clock_t start = clock();
	homework3();
	clock_t end = clock();
	cout << "�������н�������ʱ" << end - start << "ms\n";
	return 0;
}