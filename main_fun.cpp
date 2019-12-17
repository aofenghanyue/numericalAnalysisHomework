#include "homework3.h"
#include <ctime>

int main() {
	clock_t start = clock();
	homework3();
	clock_t end = clock();
	cout << "程序运行结束，用时" << end - start << "ms\n";
	return 0;
}