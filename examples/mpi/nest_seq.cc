#include <iostream>

int main() {
    double sum = 0.0;
    for (int i=0; i<10; i++) {
      for (int j=0; j<i; j++) {
	sum += j;
      }
    }
    std::cout << sum << std::endl;
    return 0;
}

